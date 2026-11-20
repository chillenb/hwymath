/*
 * AVX2/FMA3 port of the AVX-512 SVML double-precision exp, high-accuracy variant (see svml_exp_d_ha.c;
 * original svml_z0_exp_d_ha.s: Copyright (C) 2023 Intel Corporation, BSD-3-Clause; this port is under
 * the same license). Same algorithm, constants and operation order, 4 lanes. AVX-512 features that
 * AVX2 lacks are replaced by exact equivalents, so the results should be bit-identical to the
 * AVX-512 version for every input (FP exception flags raised by the vector path may differ):
 *
 *   {rz-sae} FMA     t = x*log2(e) + S is computed with the default rounding (to nearest), and
 *                    fma(x, log2(e), S - t) gives the sign of the rounding error exactly (S - t is
 *                    exact by Sterbenz); where t was rounded up, t is decremented by one ulp
 *                    (= 1/16, done on the bit pattern), which is the round-toward-zero result.
 *   vpermt2pd        16-entry lookups with vpermps + blendv (see lookup16), index j = low 4 bits of t.
 *   vscalefpd        for all lanes not handled by the scalar routine, 2^floor(N) is a normal number and
 *                    is built from the bits of t: the shifter's low bits 0x3ff0 are 1023 << 4, so
 *                    ((bits(t) - 0x42f8000000000000) >> 4) << 52 is the bit pattern of 2^M. The product
 *                    y * 2^M is then exact-or-correctly-rounded like vscalefpd.
 *   mask registers   compare + movemask; NaN inputs also go to the scalar routine (which returns
 *                    x*x, the same quieted NaN as the AVX-512 path), because the bit trick above
 *                    would not propagate NaN.
 *
 * Requires the default rounding mode (round to nearest) in MXCSR, like most SIMD math code.
 */
#include <immintrin.h>
#include <stdint.h>

#include "svml_exp_d_ha.h"

#if defined(__FAST_MATH__)
#error "svml_exp_d_ha_avx2.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif

#include "constants.inc"

#define TARGET __attribute__((target("avx2,fma")))
#define V1(c) _mm256_set1_pd(c)

/*
 * 16-entry table lookup without gathers (gathers are very slow on Skylake-derived CPUs with the
 * "Gather Data Sampling" microcode mitigation): the table sits in four registers; vpermps picks entry
 * (j & 3) of each (a double is a pair of floats, so the index pair is 2*(j&3), 2*(j&3)+1), and two
 * levels of blendv on bits 2 and 3 of j (moved into the sign bit) pick the register.
 */
TARGET static inline __m256d
lookup16(const double *tab, __m256i pair_idx, __m256d sel2, __m256d sel3)
{
    const __m256d a0 = _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab)), pair_idx));
    const __m256d a1 = _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab + 4)), pair_idx));
    const __m256d a2 = _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab + 8)), pair_idx));
    const __m256d a3 = _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab + 12)), pair_idx));
    return _mm256_blendv_pd(_mm256_blendv_pd(a0, a1, sel2), _mm256_blendv_pd(a2, a3, sel2), sel3);
}

static inline __attribute__((always_inline, target("avx2,fma"))) __m256d
exp4_ha_avx2(__m256d x)
{
    const __m256d S = V1(vShifter);
    const __m256d L = V1(vLog2e);

    /* t = RZ(x*log2(e) + S): RN first, then step back one ulp where RN rounded up */
    __m256d t = _mm256_fmadd_pd(x, L, S);
    const __m256d err = _mm256_fmadd_pd(x, L, _mm256_sub_pd(S, t));  /* exact sign of x*L + S - t */
    const __m256i up = _mm256_castpd_si256(_mm256_cmp_pd(err, _mm256_setzero_pd(), _CMP_LT_OQ));
    const __m256i tb = _mm256_add_epi64(_mm256_castpd_si256(t), up);  /* up is -1 (all ones) or 0 */
    t = _mm256_castsi256_pd(tb);

    const __m256d ax = _mm256_and_pd(x, _mm256_castsi256_pd(_mm256_set1_epi64x((long long)vAbsMaskBits)));
    const __m256d n = _mm256_sub_pd(t, S);
    /* table index j = low 4 bits of t */
    const __m256i lo2 = _mm256_slli_epi64(_mm256_and_si256(tb, _mm256_set1_epi64x(3)), 1);  /* 2*(j&3) */
    const __m256i pair_idx = _mm256_or_si256(lo2, _mm256_slli_epi64(_mm256_add_epi64(lo2, _mm256_set1_epi64x(1)), 32));
    const __m256d sel2 = _mm256_castsi256_pd(_mm256_slli_epi64(tb, 61));  /* bit 2 of j -> sign bit */
    const __m256d sel3 = _mm256_castsi256_pd(_mm256_slli_epi64(tb, 60));  /* bit 3 of j -> sign bit */
    const __m256d tlr = lookup16(vTlr, pair_idx, sel2, sel3);
    const __m256d th = lookup16(vTh, pair_idx, sel2, sel3);

    __m256d r = _mm256_fnmadd_pd(n, V1(vLn2Hi), x);  /* x - N*ln2_hi */
    /* lanes for the scalar routine: |x| >= threshold, or NaN */
    const int callout = _mm256_movemask_pd(_mm256_cmp_pd(ax, V1(vThreshold), _CMP_NLT_UQ));
    r = _mm256_fnmadd_pd(V1(vLn2Lo), n, r);          /* - N*ln2_lo */
    r = _mm256_and_pd(r, _mm256_castsi256_pd(_mm256_set1_epi64x((long long)vRMaskBits)));

    const __m256d r2 = _mm256_mul_pd(r, r);
    const __m256d p56 = _mm256_fmadd_pd(V1(vC6), r, V1(vC5));
    const __m256d p34 = _mm256_fmadd_pd(V1(vC4), r, V1(vC3));
    const __m256d p12 = _mm256_fmadd_pd(V1(vC2), r, V1(vC1));
    const __m256d p36 = _mm256_fmadd_pd(p56, r2, p34);
    const __m256d p16 = _mm256_fmadd_pd(p36, r2, p12);
    const __m256d p = _mm256_fmadd_pd(p16, r, tlr);
    const __m256d y = _mm256_fmadd_pd(th, p, th);    /* Th*(1+p) */

    /* 2^M from the bits of t (exact for every lane that does not take the callout) */
    const __m256i e = _mm256_srli_epi64(_mm256_sub_epi64(tb, _mm256_set1_epi64x(0x42f8000000000000LL)), 4);
    const __m256d scale = _mm256_castsi256_pd(_mm256_slli_epi64(e, 52));
    __m256d res = _mm256_mul_pd(y, scale);

    if (callout) {
        double xin[4], out[4];
        _mm256_storeu_pd(xin, x);
        _mm256_storeu_pd(out, res);
        for (int i = 0; i < 4; i++) {
            if (callout & (1 << i)) {
                svml_dexp_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm256_loadu_pd(out);
    }
    return res;
}

TARGET __m256d
svml_exp4_ha_avx2(__m256d x)
{
    return exp4_ha_avx2(x);
}

TARGET void
svml_exp_ha_array_avx2(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 4 <= len; i += 4) {
        _mm256_storeu_pd(y + i, exp4_ha_avx2(_mm256_loadu_pd(x + i)));
    }
    if (i < len) {
        double xin[4] = {0, 0, 0, 0}, out[4];
        for (size_t k = 0; i + k < len; k++) xin[k] = x[i + k];
        _mm256_storeu_pd(out, exp4_ha_avx2(_mm256_loadu_pd(xin)));
        for (size_t k = 0; i + k < len; k++) y[i + k] = out[k];
    }
}
