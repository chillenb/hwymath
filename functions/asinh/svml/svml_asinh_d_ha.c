/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision asinh, high-accuracy variant
 * (svml_z0_asinh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with a = |x| and asinh(-x) = -asinh(x):
 *   a < 0.0221 (kSmallX):  asinh(a) = a + a*x^2*(C1 + C2 x^2 + C3 x^4)
 *   otherwise asinh(a) = log(Xin), Xin = a + sqrt(1 + x^2) as a double-double Xin0 + Xin_low:
 *     Y = 1 + x^2 = Yh + Yl exactly (Yl from max/min(x^2, 1) and the FMA error of x^2)
 *     R0 = rsqrt14(Yh) (VRSQRT14PD, 14-bit approximation), Sh = Yh*R0, Sl = (Yh*R0)_low + R0*Yl,
 *     Eh = 1 - Sh*R0 - Sl*R0 (relative error of R0^2*Y, ~2^-14), and
 *     sqrt(Y) = (Sh + Sl)/sqrt(1 - Eh) ~ Sh + Sl + Sh*Eh*(1/2 + 3/8 Eh + 5/16 Eh^2 + 35/128 Eh^3)
 *     a >= 2^511 (kLargeX): Xin = a/8 instead (x^2 would overflow), log(2a) = log(a/8) + 4 ln2
 *   log(Xin): Rcp = rcp14(Xin0) (VRCP14PD) rounded to 1+4 mantissa bits by integer add + mask,
 *     R = Rcp*Xin0 - 1 + Rcp*Xin_low (|R| < ~2^-5), K = getexp(Rcp) (VGETEXPPD),
 *     log(Xin) = -K*ln2 + Th[k] + Tl[k] + R + R^2*P(R) with 16-entry tables -log(1 + k/16) (VPERMT2PD)
 *     and a degree-8 polynomial P, summed as a double-double with the low parts.
 * Lanes with |x| > DBL_MAX (+-inf) or NaN call the scalar routine, which returns x or x*x.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_asinh8_ha() needs AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../asinh_d.h"

#if defined(__FAST_MATH__)
#error "svml_asinh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_asinh_d_ha_constants.inc"

/* __svml_dasinh_ha_cout_rare_internal: only reached for +-inf and NaN */
int
svml_dasinh_ha_cout_rare(const double *px, double *pr)
{
    uint64_t xb;
    memcpy(&xb, px, sizeof xb);
    if (((xb >> 52) & 0x7ff) == 0x7ff && (xb & 0x000fffffffffffffULL) != 0) {  /* NaN */
        const double x = *px;
        *pr = x * x;  /* quiet NaN, raises invalid for a signaling NaN */
        return 0;
    }
    *pr = *px;  /* asinh(+-inf) = +-inf (and any other input is returned unchanged) */
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define I64(c) _mm512_set1_epi64((long long)(c))
#define AND(a, b) _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
#define XOR(a, b) _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(a), _mm512_castpd_si512(b)))
/* keeps the compiler from merging bitwise operations into vpternlogq (results are the same either way;
 * this only keeps the instruction mix identical to the assembly for make asinh-insns) */
#define OPAQUE(v) __asm__("" : "+v"(v))

/* Operations in the order of the assembly; names after its comments. */
__attribute__((target("avx512f"))) __m512d
svml_asinh8_ha(__m512d x)
{
    const __m512d one = V1(kOne);
    const __m512d x2 = _mm512_mul_round_pd(x, x, RN);
    const __m512d A = _mm512_max_round_pd(one, x2, SAE);                      /* max(x^2, 1) */
    __m512d B = _mm512_min_round_pd(one, x2, SAE);                            /* min(x^2, 1) */
    __m512d ps = _mm512_fmadd_round_pd(V1(kSmallC3), x2, V1(kSmallC2), RN);   /* small-input polynomial */
    const __m512d Y = _mm512_add_round_pd(x2, one, RN);                       /* Yh = 1 + x^2 */
    __m512d ax = AND(x, _mm512_castsi512_pd(I64(kAbsMask)));
    OPAQUE(ax);
    ps = _mm512_fmadd_round_pd(x2, ps, V1(kSmallC1), RN);
    const __m512d R0 = _mm512_rsqrt14_pd(Y);
    const __mmask8 notsmall = _mm512_cmp_round_pd_mask(ax, V1(kSmallX), _CMP_NLT_UQ, SAE);
    const __m512d Bh = _mm512_sub_round_pd(Y, A, RN);                         /* B_high */
    __m512d sign = XOR(ax, x);
    OPAQUE(sign);    const __m512d Sh = _mm512_mul_round_pd(Y, R0, RN);                        /* Sh ~ sqrt(1+x^2) */
    __m512d res = _mm512_mul_round_pd(ps, x2, RN);
    B = _mm512_sub_round_pd(B, Bh, RN);                                       /* B_low */
    const __m512d YR0l = _mm512_fmsub_round_pd(R0, Y, Sh, RN);                /* (Yh*R0)_low */
    res = _mm512_fmadd_round_pd(ax, res, ax, RN);                             /* small inputs: a + a*x^2*P */
    const __m512d x2l = _mm512_fmsub_round_pd(x, x, x2, RN);                  /* (x^2)_low */
    __m512d Sl = _mm512_add_round_pd(B, x2l, RN);                             /* Yl = (x^2)_low + B_low */
    Sl = _mm512_fmadd_round_pd(R0, Sl, YR0l, RN);                             /* Sl = (Yh*R0)_low + R0*Yl */
    const __mmask8 large = _mm512_cmp_round_pd_mask(ax, V1(kLargeX), _CMP_NLT_UQ, SAE);
    const __mmask8 callout = _mm512_cmp_round_pd_mask(ax, V1(kMaxX), _CMP_NLE_UQ, SAE);
    __m512d Eh = _mm512_fnmadd_round_pd(R0, Sh, one, RN);                     /* Eh = 1 - Sh*R0 */
    Eh = _mm512_fnmadd_round_pd(R0, Sl, Eh, RN);                              /* Eh = (1 - Sh*R0) - Sl*R0 */
    const __m512d ShEh = _mm512_mul_round_pd(Sh, Eh, RN);
    __m512d pS = _mm512_fmadd_round_pd(V1(kSqrtC4), Eh, V1(kSqrtC3), RN);     /* poly_s */
    pS = _mm512_fmadd_round_pd(Eh, pS, V1(kSqrtC2), RN);
    pS = _mm512_fmadd_round_pd(Eh, pS, V1(kSqrtC1), RN);
    const __m512d Shx = _mm512_add_round_pd(Sh, ax, RN);                      /* Sh + x */
    __m512d Sc = _mm512_fmadd_round_pd(pS, ShEh, Sl, RN);                     /* Sl + Sh*Eh*poly_s */
    const __m512d Xh = _mm512_sub_round_pd(Shx, Sh, RN);
    __m512d Xin = _mm512_add_round_pd(Shx, Sc, RN);                           /* Xin0 ~ x + sqrt(1+x^2) */
    const __m512d Xl = _mm512_sub_round_pd(ax, Xh, RN);
    const __m512d Slh = _mm512_sub_round_pd(Xin, Shx, RN);                    /* Sl_high */
    Xin = _mm512_mask_mul_round_pd(Xin, large, ax, V1(kLargeScale), RN);      /* fixup for very large inputs */
    Sc = _mm512_sub_round_pd(Sc, Slh, RN);                                    /* Sl_l */
    const __m512d rcp = _mm512_rcp14_pd(Xin);
    __m512d Xinl = _mm512_add_round_pd(Sc, Xl, RN);                           /* Xin_low */
    /* round the reciprocal to 1+4 mantissa bits */
    const __m512i rcpr = _mm512_add_epi64(_mm512_castpd_si512(rcp), I64(kRcpRoundAdd));
    Xinl = _mm512_castsi512_pd(_mm512_mask_xor_epi64(_mm512_castpd_si512(Xinl), large,
                                                     _mm512_castpd_si512(Xinl), _mm512_castpd_si512(Xinl)));
    const __m512d Rcp = AND(_mm512_castsi512_pd(rcpr), _mm512_castsi512_pd(I64(kRcpRoundMask)));
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(Rcp), 48);     /* table index */
    const __m512d Rh0 = _mm512_fmsub_round_pd(Xin, Rcp, one, RN);             /* Rcp*Xin - 1 */
    __m512d K = _mm512_getexp_round_pd(Rcp, SAE);                             /* exponents */
    const __m512d Tl = _mm512_permutex2var_pd(_mm512_loadu_pd(kLogTl), idx, _mm512_loadu_pd(kLogTl + 8));
    K = _mm512_mask_sub_round_pd(K, large, K, V1(kLargeK), RN);
    __m512d Th = _mm512_loadu_pd(kLogTh);
    const __m512d R = _mm512_fmadd_round_pd(Xinl, Rcp, Rh0, RN);              /* reduced argument */
    const __m512d Tlx = _mm512_fmadd_round_pd(Xinl, Rcp, Tl, RN);             /* Tl + Rcp*Xin_low */
    Th = _mm512_permutex2var_pd(Th, idx, _mm512_loadu_pd(kLogTh + 8));
    __m512d p = _mm512_fmadd_round_pd(V1(kLogC9), R, V1(kLogC8), RN);         /* polynomial */
    Th = _mm512_fnmadd_round_pd(V1(kL2H), K, Th, RN);                         /* -K*L2H + Th */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC7), RN);
    const __m512d Hi = _mm512_add_round_pd(Rh0, Th, RN);                      /* (K*L2H+Th) + Rh0 */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC6), RN);
    const __m512d Lo = _mm512_fnmadd_round_pd(V1(kL2L), K, Tlx, RN);          /* -K*L2L + Tl + Rcp*Xin_low */
    const __m512d nRh = _mm512_sub_round_pd(Th, Hi, RN);                      /* -R_high */
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC5), RN);
    const __m512d Rl = _mm512_add_round_pd(Rh0, nRh, RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC4), RN);
    const __m512d R2 = _mm512_mul_round_pd(R, R, RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC3), RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC2), RN);
    p = _mm512_fmadd_round_pd(R, p, V1(kLogC1), RN);
    p = _mm512_fmadd_round_pd(R2, p, Lo, RN);                                 /* (K*L2L + Tl) + R^2*Poly */
    p = _mm512_add_round_pd(p, Rl, RN);
    res = _mm512_mask_add_round_pd(res, notsmall, Hi, p, RN);
    res = XOR(res, sign);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dasinh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f"))) void
svml_asinh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_asinh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_asinh8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
