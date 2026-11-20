/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision exp, high-accuracy variant
 * (svml_z0_exp_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes):
 *   t  = x*log2(e) + 1.5*2^48       FMA rounding toward zero: N = t - 1.5*2^48 is x*log2(e)
 *                                   truncated to a multiple of 1/16, and the low 4 bits of t are j
 *   N  = M + j/16
 *   r  = x - N*ln2_hi - N*ln2_lo    (Cody-Waite reduction, 0 <= r < ln2/16 for normal inputs)
 *   2^(j/16) = Th[j]*(1 + Tlr[j])   (16-entry tables, looked up with two-source permutes)
 *   p  = r*(c1 + c2 r + ... + c6 r^5) + Tlr[j]
 *   exp(x) = scalef(Th[j] + Th[j]*p, N) = 2^floor(N) * Th[j]*(1+p)
 * Lanes with |x| >= 707.703... (and +-inf) are recomputed by a scalar routine that also signals
 * overflow/underflow; the vector result is already correct there, the callout only adds FP flags.
 * NaN compares unordered and never takes the callout.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_exp8_ha() needs AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "svml_exp_d_ha.h"

#if defined(__FAST_MATH__)
#error "svml_exp_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "constants.inc"

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
 * __svml_dexp_ha_cout_rare_internal: scalar exp for the lanes flagged by the vector path, with a
 * 64-entry table. Returns 0, or 3 (overflow) / 4 (underflow) as an error code (unused by the vector
 * routine; VML uses it for errno). Statement order mirrors the assembly.
 */
int
svml_dexp_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t e = (uint32_t)(xb >> 52) & 0x7ff;

    if (e == 0x7ff) {  /* inf or NaN */
        if ((xb >> 63) && (xb & 0x000fffffffffffffULL) == 0) {
            *pr = sZero;  /* exp(-inf) = +0 */
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
    if (x < sUnderflowX) {  /* result rounds to 0: raises underflow */
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

    /* denormal result: work with 2^(M+60), then scale by 2^-60 */
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
#define RZ (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)

__attribute__((target("avx512f"))) __m512d
svml_exp8_ha(__m512d x)
{
    /* t = x*log2(e) + 1.5*2^48 (round toward zero); its low 4 mantissa bits are the table index */
    const __m512d t = _mm512_fmadd_round_pd(x, V1(vLog2e), V1(vShifter), RZ);
    const __m512d ax = _mm512_castsi512_pd(
            _mm512_and_si512(_mm512_castpd_si512(x), _mm512_set1_epi64((long long)vAbsMaskBits)));
    const __m512d n = _mm512_sub_round_pd(t, V1(vShifter), RN);
    const __m512i idx = _mm512_castpd_si512(t);
    const __m512d tlr = _mm512_permutex2var_pd(_mm512_loadu_pd(vTlr), idx, _mm512_loadu_pd(vTlr + 8));
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));

    __m512d r = _mm512_fnmadd_round_pd(n, V1(vLn2Hi), x, RN);  /* x - N*ln2_hi */
    const __mmask8 callout = _mm512_cmp_round_pd_mask(ax, V1(vThreshold), _CMP_GE_OQ, _MM_FROUND_NO_EXC);
    r = _mm512_fnmadd_round_pd(V1(vLn2Lo), n, r, RN);            /* - N*ln2_lo */
    /* clear the top exponent bit so that |r| < 2 even for special inputs */
    r = _mm512_castsi512_pd(
            _mm512_and_si512(_mm512_castpd_si512(r), _mm512_set1_epi64((long long)vRMaskBits)));

    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __m512d p56 = _mm512_fmadd_round_pd(V1(vC6), r, V1(vC5), RN);
    const __m512d p34 = _mm512_fmadd_round_pd(V1(vC4), r, V1(vC3), RN);
    const __m512d p12 = _mm512_fmadd_round_pd(V1(vC2), r, V1(vC1), RN);
    const __m512d p36 = _mm512_fmadd_round_pd(p56, r2, p34, RN);
    const __m512d p16 = _mm512_fmadd_round_pd(p36, r2, p12, RN);
    const __m512d p = _mm512_fmadd_round_pd(p16, r, tlr, RN);    /* (exp(r)-1)(1+Tlr) + Tlr, roughly */
    const __m512d y = _mm512_fmadd_round_pd(th, p, th, RN);      /* Th*(1+p) */
    __m512d res = _mm512_scalef_round_pd(y, n, RN);              /* * 2^floor(N) = 2^M */

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dexp_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f"))) void
svml_exp_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_exp8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_exp8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
