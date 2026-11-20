/*
 * Faster AVX2/FMA3 double-precision exp, derived from the SVML high-accuracy algorithm (svml_exp_d_ha.c)
 * but no longer bit-identical to it:
 *
 *   N = M + j/T = RN(x*log2(e)) to a multiple of 1/T      one FMA with the shifter 1.5*2^(52-K), T = 2^K
 *   r = x - N*ln2_hi - N*ln2_lo, |r| <= ln2/(2T)          (the AVX-512 version truncates: 0 <= r < ln2/16)
 *   2^(j/T) = Th[j]*(1 + Tlr[j])                          T = 4: one vpermps; T = 8: two vpermps + a blend
 *   p = r*(c1 + c2 r + ...) + Tlr[j]                      minimax fit on |r| <= ln2/(2T) (tools/fit_exp_k3.py)
 *   exp(x) = (Th[j] + Th[j]*p) * 2^M                      2^M from the bits of t, as in the AVX2 port
 *
 * Differences to the bit-identical AVX2 port: half-size tables (lookup = 4 vpermps + 2 blends instead of
 * 8 + 6), no round-toward-zero emulation (saves 2 FMAs, a compare and an add), one more polynomial term.
 * Lanes with |x| >= 707.703... and NaN use the same scalar routine (svml_dexp_ha_cout_rare).
 * Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.
 */
#include <immintrin.h>
#include <stdint.h>

#include "svml_exp_d_ha.h"

#if defined(__FAST_MATH__)
#error "svml_exp_d_avx2_fast.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif

#include "constants.inc"     /* vLog2e, vLn2Hi, vLn2Lo, vThreshold, masks (shared with the HA version) */
#include "constants_fast.inc"  /* fTh, fTlr, fC1.., fShifter, fK: from tools/fit_exp_k3.py K DEG */

#define V1(c) _mm256_set1_pd(c)

/* special lanes, kept out of line so that the hot path needs no stack frame */
static __attribute__((noinline, cold, target("avx2,fma"))) __m256d
exp4_fast_callout(__m256d x, __m256d res, int mask)
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
exp4_fast(__m256d x)
{
    const __m256d S = V1(fShifter);
    const __m256d t = _mm256_fmadd_pd(x, V1(vLog2e), S);  /* RN to a multiple of 1/8; low 3 bits = j */
    const __m256i tb = _mm256_castpd_si256(t);
    const __m256d n = _mm256_sub_pd(t, S);

    /* table lookups: vpermps uses only the low 3 bits of each index, so the dword pair (2j, 2j+1) mod 8
     * is (tb << 1) duplicated into both halves of the qword, plus (0, 1) */
    const __m256i idx = _mm256_add_epi32(_mm256_shuffle_epi32(_mm256_slli_epi64(tb, 1), 0xA0),
                                         _mm256_set_epi32(1, 0, 1, 0, 1, 0, 1, 0));
#define PERM(tab) _mm256_castps_pd(_mm256_permutevar8x32_ps(_mm256_castpd_ps(_mm256_loadu_pd(tab)), idx))
#if fK == 2  /* 4 entries: one register, one vpermps */
    const __m256d th = PERM(fTh);
    const __m256d tlr = PERM(fTlr);
#elif fK == 3  /* 8 entries: two registers, bit 2 of j picks one */
    const __m256d sel = _mm256_castsi256_pd(_mm256_slli_epi64(tb, 61));
    const __m256d th = _mm256_blendv_pd(PERM(fTh), PERM(fTh + 4), sel);
    const __m256d tlr = _mm256_blendv_pd(PERM(fTlr), PERM(fTlr + 4), sel);
#else
#error "fK must be 2 or 3"
#endif
#undef PERM

    __m256d r = _mm256_fnmadd_pd(n, V1(vLn2Hi), x);
    r = _mm256_fnmadd_pd(V1(vLn2Lo), n, r);

    /* P(r) = c1 + c2 r + ... (degree 6 or 7, Estrin), p = r*P + Tlr */
    const __m256d r2 = _mm256_mul_pd(r, r);
    const __m256d r4 = _mm256_mul_pd(r2, r2);
    const __m256d q12 = _mm256_fmadd_pd(V1(fC2), r, V1(fC1));
    const __m256d q34 = _mm256_fmadd_pd(V1(fC4), r, V1(fC3));
    const __m256d q56 = _mm256_fmadd_pd(V1(fC6), r, V1(fC5));
    const __m256d a = _mm256_fmadd_pd(q34, r2, q12);
#ifdef fC8
    const __m256d b = _mm256_fmadd_pd(_mm256_fmadd_pd(V1(fC8), r, V1(fC7)), r2, q56);
#else
    const __m256d b = _mm256_fmadd_pd(V1(fC7), r2, q56);
#endif
    const __m256d P = _mm256_fmadd_pd(b, r4, a);
    const __m256d p = _mm256_fmadd_pd(P, r, tlr);
    const __m256d y = _mm256_fmadd_pd(th, p, th);

    const __m256i e = _mm256_srli_epi64(_mm256_sub_epi64(tb, _mm256_set1_epi64x(fShifterBase)), fK);
    __m256d res = _mm256_mul_pd(y, _mm256_castsi256_pd(_mm256_slli_epi64(e, 52)));

    const __m256d ax = _mm256_and_pd(x, _mm256_castsi256_pd(_mm256_set1_epi64x((long long)vAbsMaskBits)));
    const int callout = _mm256_movemask_pd(_mm256_cmp_pd(ax, V1(vThreshold), _CMP_NLT_UQ));
    if (__builtin_expect(callout, 0)) {
        res = exp4_fast_callout(x, res, callout);
    }
    return res;
}

__attribute__((target("avx2,fma"))) __m256d
svml_exp4_fast_avx2(__m256d x)
{
    return exp4_fast(x);
}

__attribute__((target("avx2,fma"))) void
svml_exp_fast_array_avx2(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 4 <= len; i += 4) {
        _mm256_storeu_pd(y + i, exp4_fast(_mm256_loadu_pd(x + i)));
    }
    if (i < len) {
        double xin[4] = {0, 0, 0, 0}, out[4];
        for (size_t k = 0; i + k < len; k++) xin[k] = x[i + k];
        _mm256_storeu_pd(out, exp4_fast(_mm256_loadu_pd(xin)));
        for (size_t k = 0; i + k < len; k++) y[i + k] = out[k];
    }
}
