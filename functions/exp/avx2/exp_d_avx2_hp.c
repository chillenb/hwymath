/*
 * AVX2/FMA3 double-precision exp aiming at Intel-HA-level accuracy (max error close to 0.5 ulp):
 * the fast variant (exp_d_avx2_fast.c) plus extra precision where its error came from.
 *
 *   N = M + j/T = RN(x*log2(e)), T = 2^K               one FMA with the shifter
 *   r_hi = x - N*Ln2Hi                                  exact: Ln2Hi has 38 significant bits, |N*T| < 2^14
 *   r_lo = -N*Ln2Lo;  r = r_hi + r_lo,  r_c = (r_hi - r) + r_lo     r + r_c = x - N*ln2 (to ~1e-25)
 *   exp(r) - 1 = r + r^2*Q(r)                           Q: minimax fit (tools/fit_poly.py K DEG hp)
 *   2^(j/T) = Th*(1 + Tlr)
 *   Th*(1+Tlr)*exp(r) = Th*(1 + r) + Th*(r^2 Q + r_c + Tlr*(1 + r)) + O(1e-33)
 *   s = fma(Th, r, Th);  e = fma(Th, r, -(s - Th))      s + e = Th + Th*r (s - Th is exact, e is the
 *                                                       rounding error of s to 2^-53 relative)
 *   y = s + fma(Th, r^2 Q + r_c + Tlr*(1+r), e)         only this last addition rounds significantly
 *   exp(x) = y * 2^M
 * Lanes with |x| >= 707.703... and NaN use the scalar routine svml_dexp_ha_cout_rare.
 * Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.
 */
#include <immintrin.h>
#include <stdint.h>

#include "../exp_d.h"

#if defined(__FAST_MATH__)
#error "exp_d_avx2_hp.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif

#include "../svml/svml_exp_d_ha_constants.inc"  /* vLog2e, vThreshold, vAbsMaskBits (shared with the HA version) */
#include "exp_d_hp_constants.inc"  /* hTh, hTlr, hC2.., hShifter, hK, hLn2Hi, hLn2Lo */

#define V1(c) _mm256_set1_pd(c)

/* special lanes, kept out of line so that the hot path needs no stack frame */
static __attribute__((noinline, cold, target("avx2,fma"))) __m256d
exp4_hp_callout(__m256d x, __m256d res, int mask)
{
    double xin[4], out[4];
    _mm256_storeu_pd(xin, x);
    _mm256_storeu_pd(out, res);
    for (int i = 0; i < 4; i++) {
        if (mask & (1 << i)) svml_dexp_ha_cout_rare(&xin[i], &out[i]);
    }
    return _mm256_loadu_pd(out);
}

static inline __attribute__((always_inline, target("avx2,fma"))) __m256d
exp4_hp(__m256d x)
{
    const __m256d S = V1(hShifter);
    const __m256d t = _mm256_fmadd_pd(x, V1(vLog2e), S);
    const __m256i tb = _mm256_castpd_si256(t);
    const __m256d n = _mm256_sub_pd(t, S);

    const __m256i idx = _mm256_add_epi32(_mm256_shuffle_epi32(_mm256_slli_epi64(tb, 1), 0xA0),
                                         _mm256_set_epi32(1, 0, 1, 0, 1, 0, 1, 0));
#define PERM(tab) _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab)), idx))
#if hK == 2
    const __m256d th = PERM(hTh);
    const __m256d tlr = PERM(hTlr);
#elif hK == 3
    const __m256d sel = _mm256_castsi256_pd(_mm256_slli_epi64(tb, 61));
    const __m256d th = _mm256_blendv_pd(PERM(hTh), PERM(hTh + 4), sel);
    const __m256d tlr = _mm256_blendv_pd(PERM(hTlr), PERM(hTlr + 4), sel);
#else
#error "hK must be 2 or 3"
#endif
#undef PERM

    /* reduction, r + r_c = x - N*ln2 */
    const __m256d rhi = _mm256_fnmadd_pd(n, V1(hLn2Hi), x);  /* exact */
    const __m256d rlo = _mm256_mul_pd(n, V1(-hLn2Lo));
    const __m256d r = _mm256_add_pd(rhi, rlo);
    const __m256d rc = _mm256_add_pd(_mm256_sub_pd(rhi, r), rlo);

    /* Q(r) = c2 + c3 r + ... (Horner in r^2 over coefficient pairs) */
    const __m256d r2 = _mm256_mul_pd(r, r);
    const __m256d q23 = _mm256_fmadd_pd(V1(hC3), r, V1(hC2));
    const __m256d q45 = _mm256_fmadd_pd(V1(hC5), r, V1(hC4));
    const __m256d q67 = _mm256_fmadd_pd(V1(hC7), r, V1(hC6));
#if defined(hC9)
    const __m256d q6 = _mm256_fmadd_pd(_mm256_fmadd_pd(V1(hC9), r, V1(hC8)), r2, q67);
#elif defined(hC8)
    const __m256d q6 = _mm256_fmadd_pd(V1(hC8), r2, q67);
#else
    const __m256d q6 = q67;
#endif
    const __m256d q = _mm256_fmadd_pd(_mm256_fmadd_pd(q6, r2, q45), r2, q23);

    /* small terms: r^2 Q + r_c + Tlr*(1 + r) */
    const __m256d small = _mm256_fmadd_pd(r2, q, _mm256_add_pd(rc, _mm256_fmadd_pd(tlr, r, tlr)));
    /* s + e = Th + Th*r */
    const __m256d s = _mm256_fmadd_pd(th, r, th);
    const __m256d e = _mm256_fmsub_pd(th, r, _mm256_sub_pd(s, th));
    const __m256d y = _mm256_add_pd(s, _mm256_fmadd_pd(th, small, e));

    const __m256i m = _mm256_srli_epi64(_mm256_sub_epi64(tb, _mm256_set1_epi64x(hShifterBase)), hK);
    __m256d res = _mm256_mul_pd(y, _mm256_castsi256_pd(_mm256_slli_epi64(m, 52)));

    const __m256d ax = _mm256_and_pd(x, _mm256_castsi256_pd(_mm256_set1_epi64x((long long)vAbsMaskBits)));
    const int callout = _mm256_movemask_pd(_mm256_cmp_pd(ax, V1(vThreshold), _CMP_NLT_UQ));
    if (__builtin_expect(callout, 0)) {
        res = exp4_hp_callout(x, res, callout);
    }
    return res;
}

__attribute__((target("avx2,fma"))) __m256d
svml_exp4_hp_avx2(__m256d x)
{
    return exp4_hp(x);
}

__attribute__((target("avx2,fma"))) void
svml_exp_hp_array_avx2(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 4 <= len; i += 4) {
        _mm256_storeu_pd(y + i, exp4_hp(_mm256_loadu_pd(x + i)));
    }
    if (i < len) {
        double xin[4] = {0, 0, 0, 0}, out[4];
        for (size_t k = 0; i + k < len; k++) xin[k] = x[i + k];
        _mm256_storeu_pd(out, exp4_hp(_mm256_loadu_pd(xin)));
        for (size_t k = 0; i + k < len; k++) y[i + k] = out[k];
    }
}
