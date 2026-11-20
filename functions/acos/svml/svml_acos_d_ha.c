/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision acos, high-accuracy variant
 * (svml_z0_acos_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), the same building blocks as asin (svml_z0_asin_d_ha.s): with a = |x|,
 *   Y = 0.5 - 0.5 a,  R = min(a^2, Y),  P: the same degree-11 polynomial, asin(s) = s + s^3 P(s^2)
 *   a < 0.5:         acos(x) = pi/2 - asin(x) = Pio2Hi + (-(x + x*R*P(R)) + Pio2Lo)
 *   a >= 0.5, x > 0: acos(x) = 2 asin(sqrt(Y)) = s + s*R*P(R),               s = 2 sqrt(Y), R = s^2/4
 *   a >= 0.5, x < 0: acos(x) = pi - 2 asin(sqrt(Y)) = PiHi + (-(s + s*R*P(R)) + PiLo)
 *     2 sqrt(Y) = Sh - Sl is the same double-double as in asin, from r = rsqrt14(Y) (vrsqrt14pd, set to 0
 *     for Y < 2^-255); the low parts (Sl, Pio2Lo/PiLo) enter before the final additions. The constants
 *     pairs are selected with masked xor/or, the sign of x is applied with xors.
 * |x| > 1 (and +-inf) are recomputed by the scalar routine (NaN result, invalid exception, return code
 * 1); NaN inputs stay on the vector path (the result is the quieted input NaN).
 *
 * The scalar routine (__svml_dacos_ha_cout_rare_internal, a complete scalar acos built like the asin
 * one) is translated mechanically, instruction by instruction, by functions/asin/tools/asm2c.py into
 * svml_dacos_ha_cout_rare.inc.
 *
 * vrsqrt14pd results are implementation-defined (Intel's are a function of the top 15 mantissa bits and
 * the exponent parity, with at most 16 significant bits); the result of the algorithm depends on them
 * (|x| >= 0.5 only). svml_acos8_ha_portable() replaces the instruction by the deterministic
 * approximation used by the Highway version, trunc16(RN(1/RN(sqrt(Y)))) (16 significant bits, relative
 * error < 2^-15.9), as a reference for the Highway targets.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only the vector functions need AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../acos_d.h"
#include "asm_model.h"

#if defined(__FAST_MATH__)
#error "svml_acos_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_acos_d_ha_constants.inc"

/* __svml_dacos_ha_cout_rare_internal */
#include "svml_dacos_ha_cout_rare.inc"

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define AVX512 __attribute__((target("avx512f")))
/* compiler barrier without instructions: keeps the instruction mix comparable with the assembly */
#define OPAQUE(v) __asm__("" : "+v"(v))

/* trunc16(RN(1/RN(sqrt(y)))): the portable replacement of vrsqrt14pd (see above) */
AVX512 static inline __m512d
rsqrt_portable(__m512d y)
{
    const __m512d r = _mm512_div_round_pd(V1(1.0), _mm512_sqrt_round_pd(y, RN), RN);
    return _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(r), _mm512_set1_epi64((long long)0xfffffff000000000ULL)));
}

AVX512 static inline __attribute__((always_inline)) __m512d
acos8(__m512d x, int portable)
{
    const __m512i signmask = _mm512_set1_epi64((long long)vSignMaskBits);
    const __m512d mx = _mm512_castsi512_pd(_mm512_or_si512(signmask, _mm512_castpd_si512(x)));  /* -|x| */
    __m512d sgn = _mm512_castsi512_pd(_mm512_and_si512(signmask, _mm512_castpd_si512(x)));
    OPAQUE(sgn);  /* keep the and and the two xors (otherwise partly fused into vpternlogq: same result) */

    /* Y = 0.5 + 0.5*(-|x|) */
    const __m512d y = _mm512_fmadd_round_pd(V1(vHalf), mx, V1(vHalf), RN);
    const __m512d x2 = _mm512_mul_round_pd(mx, mx, RN);
    __m512d r = portable ? rsqrt_portable(y) : _mm512_rsqrt14_pd(y);
    const __mmask8 tiny = _mm512_cmp_round_pd_mask(y, V1(vSmallY), _CMP_LT_OQ, SAE);
    const __mmask8 callout = _mm512_cmp_round_pd_mask(mx, V1(vMOne), _CMP_LT_OQ, SAE);  /* -|x| < -1 */
    const __m512d rr = _mm512_min_round_pd(x2, y, SAE);                                  /* R */
    r = _mm512_castsi512_pd(_mm512_mask_xor_epi64(_mm512_castpd_si512(r), tiny, _mm512_castpd_si512(r),
                                                  _mm512_castpd_si512(r)));  /* r = 0 for tiny Y */
    const __m512d y2 = _mm512_add_round_pd(y, y, RN);
    const __mmask8 big = _mm512_cmp_round_pd_mask(rr, y, _CMP_NLT_UQ, SAE);  /* R = Y: a >= 0.5 (or NaN) */
    const __mmask8 neg = _mm512_cmp_round_pd_mask(rr, x, _CMP_NLT_UQ, SAE);  /* x <= R: x < 0 if big */

    /* 2 sqrt(Y) = Sh - Sl */
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __m512d sh = _mm512_mul_round_pd(y2, r, RN);
    const __m512d e = _mm512_fmsub_round_pd(r2, y2, V1(vTwo), RN);    /* 2Y r^2 - 2 */
    const __m512d sherr = _mm512_fmsub_round_pd(r, y2, sh, RN);       /* 2Y r - Sh */
    __m512d sp = _mm512_fmadd_round_pd(V1(vSqrtC3), e, V1(vSqrtC2), RN);
    const __m512d she = _mm512_mul_round_pd(sh, e, RN);
    sp = _mm512_fmadd_round_pd(e, sp, V1(vSqrtC1), RN);
    sp = _mm512_fmadd_round_pd(e, sp, V1(vSqrtC0), RN);
    const __m512d q910 = _mm512_fmadd_round_pd(V1(vC10), rr, V1(vC9), RN);
    __m512d sl = _mm512_fmsub_round_pd(she, sp, sherr, RN);           /* Sh*e*poly(e) - Sh_err */
    __m512d q1112 = _mm512_fmadd_round_pd(V1(vC12), rr, V1(vC11), RN);
    const __m512d q34 = _mm512_fmadd_round_pd(V1(vC4), rr, V1(vC3), RN);
    const __m512d rr2 = _mm512_mul_round_pd(rr, rr, RN);
    __m512d q78 = _mm512_fmadd_round_pd(V1(vC8), rr, V1(vC7), RN);

    /* polynomial */
    q1112 = _mm512_fmadd_round_pd(rr2, q1112, q910, RN);
    const __m512d rr4 = _mm512_mul_round_pd(rr2, rr2, RN);
    const __m512d q56 = _mm512_fmadd_round_pd(V1(vC6), rr, V1(vC5), RN);
    q78 = _mm512_fmadd_round_pd(rr2, q78, q56, RN);
    __m512d q = _mm512_fmadd_round_pd(rr4, q1112, q78, RN);
    q = _mm512_fmadd_round_pd(rr2, q, q34, RN);
    sl = _mm512_mask_blend_pd(big, V1(vZero), sl);  /* Sl, or 0 */
    q = _mm512_fmadd_round_pd(rr, q, V1(vC2), RN);
    q = _mm512_fmadd_round_pd(rr, q, V1(vC1), RN);
    const __mmask8 bigneg = big & neg;
    const __m512d pr = _mm512_mul_round_pd(q, rr, RN);  /* R*P(R) */

    /* (hi, lo) = pi/2 (a < 0.5), pi (a >= 0.5, x < 0) or 0 (a >= 0.5, x > 0) */
    __m512i lo = _mm512_castpd_si512(V1(vPio2Lo));
    __m512i hi = _mm512_castpd_si512(V1(vPio2Hi));
    lo = _mm512_mask_xor_epi64(lo, big, lo, lo);
    hi = _mm512_mask_xor_epi64(hi, big, hi, hi);
    lo = _mm512_mask_or_epi64(lo, bigneg, lo, _mm512_set1_epi64((long long)vPiLoBits));
    hi = _mm512_mask_or_epi64(hi, bigneg, hi, _mm512_set1_epi64((long long)vPiHiBits));
    const __m512d slo = _mm512_castsi512_pd(_mm512_xor_si512(lo, _mm512_castpd_si512(sgn)));
    const __m512d head = _mm512_mask_blend_pd(big, mx, sh);  /* Sh, or -|x| */
    const __m512d low = _mm512_sub_round_pd(slo, sl, RN);
    const __m512d s = _mm512_sub_round_pd(head, sl, RN);     /* 2 sqrt(Y), or -|x| */
    __m512d res = _mm512_fmadd_round_pd(s, pr, low, RN);
    res = _mm512_add_round_pd(res, head, RN);
    res = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(res), _mm512_castpd_si512(sgn)));
    res = _mm512_add_round_pd(res, _mm512_castsi512_pd(hi), RN);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dacos_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

AVX512 __m512d
svml_acos8_ha(__m512d x)
{
    return acos8(x, 0);
}

AVX512 __m512d
svml_acos8_ha_portable(__m512d x)
{
    return acos8(x, 1);
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
ARRAY(svml_acos_ha_array, svml_acos8_ha)
ARRAY(svml_acos_ha_portable_array, svml_acos8_ha_portable)

#endif /* __x86_64__ */
