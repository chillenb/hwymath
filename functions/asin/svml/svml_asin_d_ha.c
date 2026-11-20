/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision asin, high-accuracy variant
 * (svml_z0_asin_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with a = |x| and the sign of x restored at the end:
 *   Y = 0.5 - 0.5 a,  R = min(a^2, Y)                          (R = a^2 for a < 0.5, Y for a >= 0.5)
 *   a < 0.5:   asin(a) = a + a*R*P(R)                           P: degree 11, asin(s) = s + s^3 P(s^2)
 *   a >= 0.5:  asin(a) = pi/2 - 2 asin(sqrt(Y)) = pi/2 - s - s*R*P(R),  s = 2 sqrt(Y), R = s^2/4
 *     2 sqrt(Y) is computed as a double-double Sh - Sl from r = rsqrt14(Y) (vrsqrt14pd, |rel err| < 2^-14,
 *     set to 0 for Y < 2^-255): Sh = RN(2Y r), Sh_err = 2Y r - Sh (exact, FMA), e = RN(2Y r^2 - 2),
 *     Sl = Sh*e*(c0 + c1 e + c2 e^2 + c3 e^3) - Sh_err   (a series for sqrt(1 + e/2)^-1 - 1, roughly);
 *     pi/2 - Sh is computed with its rounding error (Fast2Sum-like), the low parts pi/2_lo + Sl - err
 *     are added to the polynomial term, and the high part last.
 * |x| > 1 (and +-inf) are recomputed by the scalar routine (NaN result, invalid exception, return code
 * 1); NaN inputs stay on the vector path (the result is the quieted input NaN).
 *
 * The scalar routine (__svml_dasin_ha_cout_rare_internal, a complete scalar asin: x + x*x for tiny x,
 * a double-double series for |x| < 0.5, and a table-driven double-double sqrt for |x| >= 0.5) is
 * translated mechanically, instruction by instruction, by functions/asin/tools/asm2c.py into
 * svml_dasin_ha_cout_rare.inc.
 *
 * vrsqrt14pd results are implementation-defined (Intel's are a function of the top 15 mantissa bits and
 * the exponent parity, with at most 16 significant bits); the result of the algorithm depends on them
 * (|x| >= 0.5 only). Intel publishes a bit-exact reference
 * (RECIP14.c), which the Highway version emulates (functions/asin/hwy/approx14-inl.h).
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only the vector functions need AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../asin_d.h"
#include "asm_model.h"

#if defined(__FAST_MATH__)
#error "svml_asin_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_asin_d_ha_constants.inc"

/* __svml_dasin_ha_cout_rare_internal */
#include "svml_dasin_ha_cout_rare.inc"

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define AVX512 __attribute__((target("avx512f")))
/* compiler barrier without instructions: keeps the instruction mix comparable with the assembly */
#define OPAQUE(v) __asm__("" : "+v"(v))

AVX512 static inline __attribute__((always_inline)) __m512d
asin8(__m512d x)
{
    __m512d ax = _mm512_castsi512_pd(
            _mm512_and_si512(_mm512_castpd_si512(x), _mm512_set1_epi64((long long)vAbsMaskBits)));
    OPAQUE(ax);
    const __m512d x2 = _mm512_mul_round_pd(ax, ax, RN);
    const __mmask8 big = _mm512_cmp_round_pd_mask(ax, V1(vHalf), _CMP_NLT_UQ, SAE);    /* a >= 0.5 (or NaN) */
    const __mmask8 callout = _mm512_cmp_round_pd_mask(V1(vOne), ax, _CMP_LT_OQ, SAE);  /* 1 < a */
    __m512d sgn = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(ax), _mm512_castpd_si512(x)));
    OPAQUE(sgn);  /* keep x ^ |x| and the final xor (otherwise rewritten as and/vpternlogq: same result) */

    /* Y = 0.5 - 0.5 a */
    const __m512d y = _mm512_fnmadd_round_pd(V1(vHalf), ax, V1(vHalf), RN);
    __m512d r = _mm512_rsqrt14_pd(y);
    const __mmask8 tiny = _mm512_cmp_round_pd_mask(y, V1(vSmallY), _CMP_LT_OQ, SAE);
    const __m512d rr = _mm512_min_round_pd(x2, y, SAE);  /* R */
    r = _mm512_castsi512_pd(_mm512_mask_xor_epi64(_mm512_castpd_si512(r), tiny, _mm512_castpd_si512(r),
                                                  _mm512_castpd_si512(r)));  /* r = 0 for tiny Y */

    /* 2 sqrt(Y) = Sh - Sl */
    const __m512d y2 = _mm512_add_round_pd(y, y, RN);
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __m512d sh = _mm512_mul_round_pd(y2, r, RN);
    const __m512d e = _mm512_fmsub_round_pd(r2, y2, V1(vTwo), RN);    /* 2Y r^2 - 2 */
    const __m512d sherr = _mm512_fmsub_round_pd(r, y2, sh, RN);       /* 2Y r - Sh */
    __m512d sp = _mm512_fmadd_round_pd(V1(vSqrtC3), e, V1(vSqrtC2), RN);
    const __m512d she = _mm512_mul_round_pd(sh, e, RN);
    sp = _mm512_fmadd_round_pd(e, sp, V1(vSqrtC1), RN);
    const __m512d q910 = _mm512_fmadd_round_pd(V1(vC10), rr, V1(vC9), RN);
    sp = _mm512_fmadd_round_pd(e, sp, V1(vSqrtC0), RN);
    const __m512d q34 = _mm512_fmadd_round_pd(V1(vC4), rr, V1(vC3), RN);
    const __m512d sl = _mm512_fmsub_round_pd(she, sp, sherr, RN);     /* Sh*e*poly(e) - Sh_err */
    __m512d q1112 = _mm512_fmadd_round_pd(V1(vC12), rr, V1(vC11), RN);
    const __m512d rr2 = _mm512_mul_round_pd(rr, rr, RN);
    __m512d q78 = _mm512_fmadd_round_pd(V1(vC8), rr, V1(vC7), RN);

    /* -2 sqrt(Y) = Sl - Sh (a >= 0.5), or a */
    const __m512d s = _mm512_mask_sub_round_pd(ax, big, sl, sh, RN);

    /* polynomial */
    q1112 = _mm512_fmadd_round_pd(rr2, q1112, q910, RN);
    const __m512d q56 = _mm512_fmadd_round_pd(V1(vC6), rr, V1(vC5), RN);
    q78 = _mm512_fmadd_round_pd(rr2, q78, q56, RN);
    const __m512d rr4 = _mm512_mul_round_pd(rr2, rr2, RN);
    __m512d q = _mm512_fmadd_round_pd(rr4, q1112, q78, RN);
    q = _mm512_fmadd_round_pd(rr2, q, q34, RN);
    q = _mm512_fmadd_round_pd(rr, q, V1(vC2), RN);
    q = _mm512_fmadd_round_pd(rr, q, V1(vC1), RN);
    __m512d pr = _mm512_mul_round_pd(q, rr, RN);   /* R*P(R) */

    /* pi/2 - Sh with its rounding error */
    const __m512d hi = _mm512_sub_round_pd(V1(vPio2Hi), sh, RN);
    const __m512d lo = _mm512_add_round_pd(V1(vPio2Lo), sl, RN);
    const __m512d t = _mm512_sub_round_pd(V1(vPio2Hi), hi, RN);
    const __m512d err = _mm512_sub_round_pd(sh, t, RN);
    /* low part: Pio2Lo + Sl - err (a >= 0.5), or 0 */
    const __m512d low = _mm512_mask_sub_round_pd(V1(vZero), big, lo, err, RN);
    pr = _mm512_fmadd_round_pd(s, pr, low, RN);
    const __m512d head = _mm512_mask_blend_pd(big, s, hi);  /* a, or pi/2 - Sh */
    __m512d res = _mm512_add_round_pd(head, pr, RN);
    res = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(res), _mm512_castpd_si512(sgn)));

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dasin_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

AVX512 __m512d
svml_asin8_ha(__m512d x)
{
    return asin8(x);
}

#define ARRAY(NAME, KERNEL)                                                                         \
    AVX512 void NAME(const double *x, double *y, size_t len)                                        \
    {                                                                                               \
        size_t i = 0;                                                                               \
        for (; i + 8 <= len; i += 8) {                                                              \
            _mm512_storeu_pd(y + i, KERNEL(_mm512_loadu_pd(x + i)));                                \
        }                                                                                           \
        if (i < len) {                                                                              \
            const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);                                   \
            _mm512_mask_storeu_pd(y + i, m, KERNEL(_mm512_maskz_loadu_pd(m, x + i)));               \
        }                                                                                           \
    }
ARRAY(svml_asin_ha_array, svml_asin8_ha)

#endif /* __x86_64__ */
