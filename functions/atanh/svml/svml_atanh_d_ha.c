/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision atanh, high-accuracy variant
 * (svml_z0_atanh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with y = |x| and atanh(-x) = -atanh(x):
 *   atanh(y) = (log(1 + y) - log(1 - y))/2, Yp = 1 + y and Ym = 1 - y as double-doubles
 *   (Yp_low = y - (Yp - 1), Ym_low = y + (Ym - 1), i.e. Ym = Ym_high - Ym_low exactly);
 *   each log as in asinh/acosh: Rcp = rcp14(Yp or Ym) (VRCP14PD) rounded to 1+4 mantissa bits,
 *   R = Rcp*Y - 1 + Rcp*Y_low, K = getexp(Rcp) (VGETEXPPD), log(1 + k/16) from 16-entry tables
 *   (VPERMT2PD/VPERMI2PD), degree-10 polynomial; the two logs are subtracted as double-doubles
 *   ((Km - Kp)*ln2 + table differences + Rp - Rm + polynomial terms) and the sum is multiplied by +-1/2.
 * Lanes with |x| >= 1 or NaN call the scalar routine.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_atanh8_ha() needs AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#include <math.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../atanh_d.h"

#if defined(__FAST_MATH__)
#error "svml_atanh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_atanh_d_ha_constants.inc"

/* __svml_datanh_ha_cout_rare_internal; the vector routine calls it for |x| >= 1 and NaN */
int
svml_datanh_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    uint64_t xb;
    memcpy(&xb, px, sizeof xb);
    if (((xb >> 52) & 0x7ff) == 0x7ff && (xb & 0x000fffffffffffffULL) != 0) {  /* NaN */
        *pr = x * x;
        return 0;
    }
    if ((uint32_t)xb == 0 && (uint32_t)(xb >> 32 & 0x7fffffff) == 0x3ff00000) {  /* x = +-1 */
        volatile double zero = kCoutTab[1];
        *pr = x / zero;  /* +-inf, raises divide-by-zero */
        return 2;
    }
    /* |x| > 1: 0*inf at run time raises invalid; the result is the x86 default NaN (negative quiet NaN),
     * returned as such on every architecture */
    volatile double zero = kCoutTab[1];
    volatile double nan = zero * kCoutTab[2];
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
#define CI(a) _mm512_castpd_si512(a)
#define CD(a) _mm512_castsi512_pd(a)
/* keeps the compiler from merging bitwise operations into vpternlogq (results are the same either way;
 * this only keeps the instruction mix identical to the assembly for make atanh-insns) */
#define OPAQUE(v) __asm__("" : "+v"(v))

/* Operations in the order of the assembly; names after its comments. */
__attribute__((target("avx512f"))) __m512d
svml_atanh8_ha(__m512d x)
{
    const __m512d one = V1(kOne);
    __m512d y = CD(_mm512_and_si512(CI(x), I64(kAbsMask)));
    OPAQUE(y);
    const __m512d Yp = _mm512_add_round_pd(y, one, RN);                       /* 1 + y */
    const __m512d Ym = _mm512_sub_round_pd(one, y, RN);                       /* 1 - y */
    __m512d sign = CD(_mm512_xor_si512(CI(x), CI(y)));
    OPAQUE(sign);
    const __m512d Yph = _mm512_sub_round_pd(Yp, one, RN);                     /* Yp_high */
    const __m512d nYmh = _mm512_sub_round_pd(Ym, one, RN);                    /* -Ym_high */
    const __m512d rcpP = _mm512_rcp14_pd(Yp);                                 /* RcpP ~ 1/Yp */
    const __m512d rcpM = _mm512_rcp14_pd(Ym);                                 /* RcpM ~ 1/Ym */
    const __mmask8 callout = _mm512_cmp_round_pd_mask(y, one, _CMP_NLT_UQ, SAE);
    /* round the reciprocals to 1+4 mantissa bits */
    const __m512i rP = _mm512_add_epi64(CI(rcpP), I64(kRcpRoundAdd));
    const __m512i rM = _mm512_add_epi64(CI(rcpM), I64(kRcpRoundAdd));
    const __m512d Ypl = _mm512_sub_round_pd(y, Yph, RN);                      /* Yp_low */
    const __m512d RcpP = CD(_mm512_and_si512(rP, I64(kRcpRoundMask)));
    const __m512d RcpM = CD(_mm512_and_si512(rM, I64(kRcpRoundMask)));
    const __m512d Yml = _mm512_add_round_pd(y, nYmh, RN);                     /* Ym_low */
    __m512d Rp = _mm512_fmsub_round_pd(RcpP, Yp, one, RN);                    /* RcpP*Yp - 1 */
    const __m512d Kp = _mm512_getexp_round_pd(RcpP, SAE);
    Rp = _mm512_fmadd_round_pd(Ypl, RcpP, Rp, RN);                            /* + RcpP*Yp_low */
    const __m512i idxM = _mm512_srli_epi64(CI(RcpM), 48);
    const __m512i idxP = _mm512_srli_epi64(CI(RcpP), 48);
    __m512d Rm = _mm512_fmsub_round_pd(Ym, RcpM, one, RN);                    /* RcpM*Ym - 1 */
    Rm = _mm512_fnmadd_round_pd(Yml, RcpM, Rm, RN);                           /* - RcpM*Ym_low */
    const __m512d Km = _mm512_getexp_round_pd(RcpM, SAE);
    const __m512d dK = _mm512_sub_round_pd(Km, Kp, RN);                       /* Km - Kp */
    const __m512d Tl0 = _mm512_loadu_pd(kLogTl), Tl1 = _mm512_loadu_pd(kLogTl + 8);
    const __m512d Th0 = _mm512_loadu_pd(kLogTh), Th1 = _mm512_loadu_pd(kLogTh + 8);
    const __m512d TlP = _mm512_permutex2var_pd(Tl0, idxP, Tl1);
    const __m512d TlM = _mm512_permutex2var_pd(Tl0, idxM, Tl1);
    const __m512d ThP = _mm512_permutex2var_pd(Th0, idxP, Th1);
    const __m512d ThM = _mm512_permutex2var_pd(Th0, idxM, Th1);
    const __m512d dTl = _mm512_sub_round_pd(TlM, TlP, RN);
    __m512d pP = _mm512_fmadd_round_pd(V1(kLogC9), Rp, V1(kLogC8), RN);       /* polynomials */
    __m512d pM = _mm512_fmadd_round_pd(V1(kLogC9), Rm, V1(kLogC8), RN);
    const __m512d dTh = _mm512_sub_round_pd(ThM, ThP, RN);                    /* table values */
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC7), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC7), RN);
    const __m512d Hi0 = _mm512_fmadd_round_pd(V1(kL2H), dK, dTh, RN);        /* K*L2H + Th */
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC6), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC6), RN);
    const __m512d Lo0 = _mm512_fmadd_round_pd(V1(kL2L), dK, dTl, RN);        /* K*L2L + Tl */
    const __m512d Rp2 = _mm512_mul_round_pd(Rp, Rp, RN);
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC5), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC5), RN);
    const __m512d H1 = _mm512_add_round_pd(Rp, Hi0, RN);                      /* (K*L2H+Th) + Rph */
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC4), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC4), RN);
    const __m512d Rm2 = _mm512_mul_round_pd(Rm, Rm, RN);
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC3), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC3), RN);
    const __m512d nRph = _mm512_sub_round_pd(Hi0, H1, RN);                    /* -Rp_high */
    const __m512d H2 = _mm512_sub_round_pd(H1, Rm, RN);                       /* (K*L2H+Th) + (Rph - Rmh) */
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC2), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC2), RN);
    pP = _mm512_fmadd_round_pd(Rp, pP, V1(kLogC1), RN);
    pM = _mm512_fmadd_round_pd(Rm, pM, V1(kLogC1), RN);
    const __m512d Rpl = _mm512_add_round_pd(Rp, nRph, RN);
    const __m512d nRmh = _mm512_sub_round_pd(H2, H1, RN);                     /* -Rm_high */
    pP = _mm512_fmadd_round_pd(Rp2, pP, Lo0, RN);                             /* (K*L2L + Tl) + Rp^2*PolyP */
    pM = _mm512_fnmadd_round_pd(Rm2, pM, Rpl, RN);                            /* Rpl - Rm^2*PolyM */
    const __m512d Rml = _mm512_add_round_pd(Rm, nRmh, RN);
    const __m512d half = CD(_mm512_or_si512(CI(sign), CI(V1(kHalf))));        /* +-1/2 */
    __m512d lo = _mm512_add_round_pd(pP, pM, RN);                             /* low part */
    lo = _mm512_sub_round_pd(lo, Rml, RN);
    __m512d res = _mm512_add_round_pd(H2, lo, RN);
    res = _mm512_mul_round_pd(res, half, RN);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_datanh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f"))) void
svml_atanh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_atanh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_atanh8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
