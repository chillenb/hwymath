/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision acosh, high-accuracy variant
 * (svml_z0_acosh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), the same as asinh's (functions/asinh) with Y = x^2 - 1:
 *   Y = x^2 - 1 = Yh + Yl exactly (from the FMA error of x^2), R0 = rsqrt14(Yh) (VRSQRT14PD),
 *   S = sqrt(Y) ~ Sh + Sl + Sh*Eh*poly_s(Eh), with Sh = Yh*R0 and Eh = 1 - R0^2*Y the relative error
 *   Yh < 2^-11 (kSmallX):  acosh(x) = asinh(S) = S + S*Y*(C1 + C2 Y + C3 Y^2)
 *   otherwise acosh(x) = log(x + S) (log(x/8) + 4 ln2 for x >= 2^511, where x^2 overflows), with the
 *   log of a double-double computed as in asinh: VRCP14PD reciprocal rounded to 1+4 mantissa bits,
 *   16-entry tables (VPERMT2PD), VGETEXPPD exponent, degree-10 polynomial.
 * Lanes with x <= 1, x > DBL_MAX (+inf) or NaN call the scalar routine (x = 1 too: there Yh = 0 and
 * rsqrt14(0) = +inf would make the vector result a NaN).
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_acosh8_ha() needs AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#include <math.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../acosh_d.h"

#if defined(__FAST_MATH__)
#error "svml_acosh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_acosh_d_ha_constants.inc"

/* __svml_dacosh_ha_cout_rare_internal; the vector routine calls it for x <= 1, x > DBL_MAX and NaN */
int
svml_dacosh_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    uint64_t xb;
    memcpy(&xb, px, sizeof xb);
    if (((xb >> 52) & 0x7ff) == 0x7ff && (xb & 0x000fffffffffffffULL) != 0) {  /* NaN */
        *pr = x * x;
        return 0;
    }
    const uint32_t lo = (uint32_t)xb, hi = (uint32_t)(xb >> 32);
    if (lo == 0) {
        if (hi == 0x3ff00000) {  /* x = 1 */
            *pr = kCoutTab[1];   /* +0 */
            return 0;
        }
        if (hi == 0x7ff00000) {  /* +inf */
            *pr = x;
            return 0;
        }
    }
    /* x < 1 (or -inf): inf*0 at run time raises invalid; the result is the x86 default NaN (negative
     * quiet NaN), returned as such on every architecture */
    volatile double inf = kCoutTab[2];
    volatile double nan = inf * kCoutTab[1];
    (void)nan;
    const uint64_t qnan = 0xfff8000000000000ULL;
    memcpy(pr, &qnan, sizeof qnan);
    return 1;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define I64(c) _mm512_set1_epi64((long long)(c))
#define AND(a, b) _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))

/* Operations in the order of the assembly; names after its comments. */
__attribute__((target("avx512f"))) __m512d
svml_acosh8_ha(__m512d x)
{
    const __m512d one = V1(kOne);
    const __m512d x2 = _mm512_mul_round_pd(x, x, RN);
    const __mmask8 le1 = _mm512_cmp_round_pd_mask(x, one, _CMP_LE_OQ, SAE);
    const __m512d Y0 = _mm512_sub_round_pd(x2, one, RN);                      /* x^2 - 1 */
    const __m512d t = _mm512_sub_round_pd(x2, Y0, RN);
    const __m512d x2l = _mm512_fmsub_round_pd(x, x, x2, RN);                  /* (x^2)_low */
    const __m512d nOnel = _mm512_sub_round_pd(t, one, RN);                    /* -One_l */
    const __m512d Yl0 = _mm512_add_round_pd(x2l, nOnel, RN);
    const __m512d Y = _mm512_add_round_pd(Y0, Yl0, RN);                       /* Yh */
    const __mmask8 large = _mm512_cmp_round_pd_mask(x, V1(kLargeX), _CMP_NLT_UQ, SAE);
    const __m512d R0 = _mm512_rsqrt14_pd(Y);
    const __mmask8 notsmall = _mm512_cmp_round_pd_mask(Y, V1(kSmallX), _CMP_NLT_UQ, SAE);
    __m512d ps = _mm512_fmadd_round_pd(V1(kSmallC3), Y, V1(kSmallC2), RN);    /* small-input polynomial */
    const __m512d Sh = _mm512_mul_round_pd(Y, R0, RN);                        /* Sh ~ sqrt(x^2 - 1) */
    ps = _mm512_fmadd_round_pd(Y, ps, V1(kSmallC1), RN);
    const __m512d Ylh = _mm512_sub_round_pd(Y, Y0, RN);
    const __m512d Shx = _mm512_add_round_pd(Sh, x, RN);                       /* Sh + x */
    const __m512d pY = _mm512_mul_round_pd(ps, Y, RN);
    const __mmask8 big = _mm512_cmp_round_pd_mask(x, V1(kMaxX), _CMP_NLE_UQ, SAE);
    const __m512d Yl = _mm512_sub_round_pd(Yl0, Ylh, RN);
    const __mmask8 callout = big | le1;
    const __m512d YR0l = _mm512_fmsub_round_pd(R0, Y, Sh, RN);                /* (Yh*R0)_low */
    __m512d Eh = _mm512_fnmadd_round_pd(R0, Sh, one, RN);                     /* Eh = 1 - Sh*R0 */
    const __m512d Sl = _mm512_fmadd_round_pd(R0, Yl, YR0l, RN);               /* Sl = (Yh*R0)_low + R0*Yl */
    Eh = _mm512_fnmadd_round_pd(R0, Sl, Eh, RN);                              /* Eh = (1 - Sh*R0) - Sl*R0 */
    const __m512d ShEh = _mm512_mul_round_pd(Sh, Eh, RN);
    __m512d pS = _mm512_fmadd_round_pd(V1(kSqrtC4), Eh, V1(kSqrtC3), RN);     /* poly_s */
    pS = _mm512_fmadd_round_pd(Eh, pS, V1(kSqrtC2), RN);
    const __m512d Shh = _mm512_sub_round_pd(Shx, x, RN);
    pS = _mm512_fmadd_round_pd(Eh, pS, V1(kSqrtC1), RN);
    const __m512d Shl = _mm512_sub_round_pd(Sh, Shh, RN);
    const __m512d Sc = _mm512_fmadd_round_pd(pS, ShEh, Sl, RN);               /* Sl + Sh*Eh*poly_s */
    __m512d s = _mm512_add_round_pd(Sh, Sc, RN);                              /* small inputs */
    __m512d Xin = _mm512_add_round_pd(Shx, Sc, RN);                           /* Xin0 ~ x + sqrt(x^2 - 1) */
    s = _mm512_fmadd_round_pd(pY, s, Sc, RN);
    const __m512d Slh = _mm512_sub_round_pd(Xin, Shx, RN);                    /* Sl_high */
    __m512d res = _mm512_add_round_pd(s, Sh, RN);                             /* small inputs: S + S*Y*P */
    const __m512d Sll = _mm512_sub_round_pd(Sc, Slh, RN);                     /* Sl_l */
    Xin = _mm512_mask_mul_round_pd(Xin, large, x, V1(kLargeScale), RN);       /* fixup for very large inputs */
    __m512d Xinl = _mm512_add_round_pd(Sll, Shl, RN);                         /* Xin_low */
    const __m512d rcp = _mm512_rcp14_pd(Xin);
    Xinl = _mm512_castsi512_pd(_mm512_mask_xor_epi64(_mm512_castpd_si512(Xinl), large,
                                                     _mm512_castpd_si512(Xinl), _mm512_castpd_si512(Xinl)));
    /* round the reciprocal to 1+4 mantissa bits */
    const __m512i rcpr = _mm512_add_epi64(_mm512_castpd_si512(rcp), I64(kRcpRoundAdd));
    const __m512d Rcp = AND(_mm512_castsi512_pd(rcpr), _mm512_castsi512_pd(I64(kRcpRoundMask)));
    const __m512d Rh0 = _mm512_fmsub_round_pd(Xin, Rcp, one, RN);             /* Rcp*Xin - 1 */
    __m512d K = _mm512_getexp_round_pd(Rcp, SAE);                             /* exponents */
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(Rcp), 48);     /* table index */
    K = _mm512_mask_sub_round_pd(K, large, K, V1(kLargeK), RN);
    const __m512d Tl = _mm512_permutex2var_pd(_mm512_loadu_pd(kLogTl), idx, _mm512_loadu_pd(kLogTl + 8));
    __m512d Th = _mm512_permutex2var_pd(_mm512_loadu_pd(kLogTh), idx, _mm512_loadu_pd(kLogTh + 8));
    const __m512d R = _mm512_fmadd_round_pd(Xinl, Rcp, Rh0, RN);              /* reduced argument */
    const __m512d Tlx = _mm512_fmadd_round_pd(Xinl, Rcp, Tl, RN);             /* Tl + Rcp*Xin_low */
    Th = _mm512_fnmadd_round_pd(V1(kL2H), K, Th, RN);                         /* -K*L2H + Th */
    __m512d p = _mm512_fmadd_round_pd(V1(kLogC9), R, V1(kLogC8), RN);         /* polynomial */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC7), RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC6), RN);
    const __m512d Lo = _mm512_fnmadd_round_pd(V1(kL2L), K, Tlx, RN);          /* -K*L2L + Tl + Rcp*Xin_low */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC5), RN);
    const __m512d Hi = _mm512_add_round_pd(Rh0, Th, RN);                      /* (K*L2H+Th) + Rh0 */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC4), RN);
    const __m512d nRh = _mm512_sub_round_pd(Th, Hi, RN);                      /* -R_high */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC3), RN);
    const __m512d R2 = _mm512_mul_round_pd(R, R, RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC2), RN);
    const __m512d Rl = _mm512_add_round_pd(Rh0, nRh, RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC1), RN);
    p = _mm512_fmadd_round_pd(R2, p, Lo, RN);                                 /* (K*L2L + Tl) + R^2*Poly */
    p = _mm512_add_round_pd(p, Rl, RN);
    res = _mm512_mask_add_round_pd(res, notsmall, Hi, p, RN);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dacosh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f"))) void
svml_acosh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_acosh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_acosh8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
