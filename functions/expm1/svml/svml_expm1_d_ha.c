/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision expm1, high-accuracy variant
 * (svml_z0_expm1_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), all operations rounded to nearest:
 *   t  = x*log2(e) + 1.5*2^48 (+1023)  N = t - shifter is x*log2(e) rounded to a multiple of 1/16, and
 *                                      the low 4 bits of t are j
 *   R  = x - N*ln2_hi (one FMA), R' = R - N*ln2_lo, Rl = (R - R') - N*ln2_lo (its rounding error)
 *   Th' = scalef(Th[j], max(N, -128)) = 2^(N clamped)*Th[j], rounded;  2^(j/16) = Th[j]*(1 + Tlr[j])
 *   expm1(x) = Th'*(1 + R + R^2*P(R) + Tlr*(1 + R)) - 1, with P(R) = c2 + c3 R + ... + c8 R^6,
 *     evaluated as (Th' - 1) + Th'*R in double-double (Fast2Sum of Th' and -1 ordered by magnitude,
 *     exact product error with an FMA), plus Th'*(R^2 P(R) + Tlr*R + Tlr + Rl) and the low parts;
 *   the sign of x is ORed into the result (so expm1(-0) = -0).
 * R, R' and Rl have their top exponent bit cleared ("|R| < 2 even for special cases": -inf, huge |x|).
 * Lanes with x >= 707.703... (and +inf, NaN) are recomputed by a scalar routine (exp(x) - 1 = exp(x) in
 * that range, with overflow handling and FP flags); large negative x needs no callout (result -1).
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_expm18_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../expm1_d.h"

#if defined(__FAST_MATH__)
#error "svml_expm1_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_expm1_d_ha_constants.inc"
#define sMinusOne (-0x1p+0)  /* an immediate (0xbff0000000000000) in the assembly */

static inline uint64_t
dbits(double d)
{
    uint64_t u;
    memcpy(&u, &d, sizeof u);
    return u;
}

static inline double
bitsd(uint64_t u)
{
    double d;
    memcpy(&d, &u, sizeof d);
    return d;
}

/*
 * __svml_dexpm1_ha_cout_rare_internal: scalar routine for the lanes flagged by the vector path. It is
 * SVML's scalar exp (as in exp's callout, same 64-entry table) preceded by x < 0 -> -1, so it is only
 * correct where exp(x) - 1 rounds to exp(x) or -1. Returns 0, or 3 (overflow) / 4 (underflow).
 * Statement order mirrors the assembly.
 */
int
svml_dexpm1_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t e = (uint32_t)(xb >> 52) & 0x7ff;

    if (x < sZero) {  /* signalling compare (comisd): raises invalid for any NaN */
        *pr = sMinusOne;
        return 0;
    }
    if (e == 0x7ff) {  /* +inf or NaN (-inf took the branch above) */
        if ((xb >> 63) && (xb & 0x000fffffffffffffULL) == 0) {
            *pr = sZero;  /* exp(-inf) = +0 (unreachable here) */
        }
        else {
            *pr = x * x;  /* +inf -> +inf, NaN -> quiet NaN */
        }
        return 0;
    }
    if (e <= 970) {  /* |x| < 2^-52 */
        *pr = sOne + x;
        return 0;
    }
    if (sOverflowX < x) {  /* x > log(DBL_MAX): +inf, raises overflow */
        volatile double h = sHuge;  /* computed at run time to raise the overflow flag */
        *pr = h * h;
        return 3;
    }
    if (x < sUnderflowX) {  /* result rounds to 0: raises underflow (unreachable here) */
        volatile double t = sTiny;  /* computed at run time to raise the underflow flag */
        *pr = t * t;
        return 4;
    }

    /* x = (M*64 + j) * ln2/64 + r */
    double t = x * s64OverLn2;
    double s = t + sShifter;
    const uint32_t ni = (uint32_t)dbits(s);  /* low 32 bits of s: N in two's complement */
    const uint32_t j = ni & 63;
    double n = s - sShifter;
    double r = x - n * sLn2o64Hi;
    r = r - n * sLn2o64Lo;

    double p = sC6 * r;
    p = p + sC5;
    p = p * r;
    const double thi = sTab[2 * j];
    p = p + sC4;
    p = p * r;
    p = p + sC3;
    p = p * r;
    p = p + sC2;
    p = p * r;
    p = p * r;
    p = p + r;
    p = p + sTab[2 * j + 1];
    p = p * thi;  /* p = Thi*(exp(r)-1 + Tlo) */

    if (!(x < sDenormX)) {  /* normal result */
        uint32_t be = ((ni >> 6) + 1023) & 2047;  /* biased exponent of 2^M */
        p = p + thi;
        if (be <= 2046) {
            *pr = p * bitsd((uint64_t)be << 52);
        }
        else {  /* M = 1024: scale by 2^1023, then by 2 */
            be = (be - 1) & 2047;
            p = p * bitsd((uint64_t)be << 52);
            *pr = p * sTwo;
        }
        return 0;
    }

    /* denormal result (unreachable here): work with 2^(M+60), then scale by 2^-60 */
    const uint32_t be = ((ni >> 6) + 1083) & 2047;
    const double scale = bitsd((uint64_t)be << 52);
    const double lo = p * scale;
    const double hi = scale * thi;
    const double a = hi + lo;
    if (be <= 50) {
        *pr = a * sTwoM60;
    }
    else {  /* split a so that the final scaling rounds only once */
        double b = hi - a;
        b = b + lo;
        const double c = a * sSplitter;
        const double d = a + c;
        double ah = d - c;
        double al = a - ah;
        al = b + al;
        ah = ah * sTwoM60;
        al = al * sTwoM60;
        *pr = ah + al;
    }
    return 4;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define BITS(c) _mm512_castsi512_pd(_mm512_set1_epi64((long long)(c)))

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_expm18_ha(__m512d x)
{
    /* t = x*log2(e) + shifter: N = t - shifter is x*log2(e) rounded to a multiple of 1/16 */
    const __m512d t = _mm512_fmadd_round_pd(x, V1(vLog2e), V1(vShifter), RN);
    const __mmask8 callout = _mm512_cmp_round_pd_mask(x, V1(vThreshold), _CMP_NLT_UQ, _MM_FROUND_NO_EXC);
    const __m512d n = _mm512_sub_round_pd(t, V1(vShifter), RN);
    const __m512i idx = _mm512_castpd_si512(t);
    __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));
    const __m512d tlr = _mm512_permutex2var_pd(_mm512_loadu_pd(vTlr), idx, _mm512_loadu_pd(vTlr + 8));
    __m512d sgn = _mm512_and_pd(x, BITS(vSignMaskBits));
    /* opaque to the optimizer: otherwise GCC fuses this AND with the final OR into one vpternlogq (same
     * result), and the instruction mix no longer matches the assembly's (make expm1-insns) */
    __asm__("" : "+v"(sgn));
    const __m512d r = _mm512_fnmadd_round_pd(n, V1(vLn2Hi), x, RN);          /* R = x - N*ln2_hi */
    const __m512d nc = _mm512_max_round_pd(n, V1(vMinN), _MM_FROUND_NO_EXC);  /* max(N, -128) */
    th = _mm512_scalef_round_pd(th, nc, RN);                                  /* Th' = Th*2^floor(N) */
    __m512d rr = _mm512_fnmadd_round_pd(n, V1(vLn2Lo), r, RN);               /* R' = R - N*ln2_lo */
    const __mmask8 small = _mm512_cmp_round_pd_mask(th, V1(vOne), _CMP_LT_OQ, _MM_FROUND_NO_EXC);
    __m512d rl = _mm512_sub_round_pd(r, rr, RN);                              /* (N*ln2_lo)_high */
    rr = _mm512_and_pd(rr, BITS(vRMaskBits));                                 /* |R'| < 2 */
    rl = _mm512_fnmadd_round_pd(V1(vLn2Lo), n, rl, RN);                       /* Rl = -(N*ln2_lo)_low */
    const __m512d p56 = _mm512_fmadd_round_pd(V1(vC6), rr, V1(vC5), RN);
    const __m512d p34 = _mm512_fmadd_round_pd(V1(vC4), rr, V1(vC3), RN);
    rl = _mm512_and_pd(rl, BITS(vRMaskBits));
    const __m512d r2 = _mm512_mul_round_pd(rr, rr, RN);
    const __m512d rltlr = _mm512_add_round_pd(rl, tlr, RN);                   /* Rl + Tlr */
    __m512d p = _mm512_fmadd_round_pd(V1(vC8), rr, V1(vC7), RN);
    const __m512d lo1 = _mm512_fmadd_round_pd(rr, tlr, rltlr, RN);            /* Rl + Tlr + Tlr*R */
    p = _mm512_fmadd_round_pd(r2, p, p56, RN);
    /* Fast2Sum of Th' and -1: the larger magnitude first */
    const __m512d big = _mm512_mask_blend_pd(small, th, V1(vMinusOne));
    const __m512d sml = _mm512_mask_blend_pd(small, V1(vMinusOne), th);
    p = _mm512_fmadd_round_pd(r2, p, p34, RN);
    const __m512d tm1 = _mm512_add_round_pd(th, V1(vMinusOne), RN);          /* Th' - 1 */
    p = _mm512_fmadd_round_pd(rr, p, V1(vC2), RN);
    const __m512d tl1 = _mm512_add_round_pd(_mm512_sub_round_pd(big, tm1, RN), sml, RN);  /* its error */
    p = _mm512_fmadd_round_pd(r2, p, lo1, RN);                                /* R^2 P(R) + Tlr*(1+R) + Rl */
    const __m512d hi = _mm512_fmadd_round_pd(rr, th, tm1, RN);                /* Th' - 1 + (Th'*R)_high */
    __m512d lo = _mm512_sub_round_pd(hi, tm1, RN);
    lo = _mm512_fmsub_round_pd(rr, th, lo, RN);                               /* (Th'*R)_low */
    __m512d res = _mm512_fmadd_round_pd(p, th, tl1, RN);                      /* Tl1 + Th'*(...) */
    res = _mm512_add_round_pd(res, lo, RN);
    res = _mm512_add_round_pd(hi, res, RN);
    res = _mm512_or_pd(res, sgn);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dexpm1_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_expm1_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_expm18_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_expm18_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
