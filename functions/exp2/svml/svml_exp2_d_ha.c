/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision exp2, high-accuracy variant
 * (svml_z0_exp2_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes):
 *   r  = vreducepd(x, imm 0x41) = x - floor(16 x)/16      in [0, 1/16); the subtraction rounds
 *                                                        toward -inf (inexact only for -1/16 < x < 0)
 *   t  = x + 1.5*2^48 rounded toward -inf               low 4 mantissa bits of t: j = floor(16 x) mod 16
 *   2^(j/16) = Th[j] + Tl[j]                            (16-entry tables, two-source permutes)
 *   y  = Th + (Tl + (Th*r) * (c1 + c2 r + ... + c6 r^5))
 *   exp2(x) = scalef(y, x) = y * 2^floor(x)
 * Lanes with |x| >= 1021.99951171875 (and +-inf) are recomputed by a scalar routine with its own
 * 64-entry table, which handles overflow, denormal results and underflow and raises the FP flags.
 * NaN compares unordered and never takes the callout: the vector path returns it quieted.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_exp28_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../exp2_d.h"

#if defined(__FAST_MATH__)
#error "svml_exp2_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_exp2_d_ha_constants.inc"

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
 * __svml_dexp2_ha_cout_rare_internal: scalar exp2 for the lanes flagged by the vector path, with a
 * 64-entry table: x = M + j/64 + r, |r| <= 1/128. Returns 0, or 3 (overflow) / 4 (underflow) as an
 * error code (unused by the vector routine). Statement order mirrors the assembly.
 */
int
svml_dexp2_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t e = (uint32_t)(xb >> 52) & 0x7ff;

    if (e == 0x7ff) {  /* inf or NaN */
        if ((xb >> 63) && (xb & 0x000fffffffffffffULL) == 0) {
            *pr = sZero;  /* exp2(-inf) = +0 */
        }
        else {
            *pr = x * x;  /* +inf -> +inf, NaN -> quiet NaN */
        }
        return 0;
    }
    if (!(x < sOverflowX)) {  /* x >= 1024: +inf, raises overflow */
        volatile double h = sHuge;  /* computed at run time to raise the overflow flag */
        *pr = h * h;
        return 3;
    }
    if (!(x > sUnderflowX)) {  /* x <= -1075: result rounds to 0, raises underflow */
        volatile double t = sTiny;  /* computed at run time to raise the underflow flag */
        *pr = t * t;
        return 4;
    }

    /* x = (M*64 + j)/64 + r */
    double s = x + sShifter;
    const uint32_t ni = (uint32_t)dbits(s);  /* low 32 bits of s: M*64 + j in two's complement */
    const uint32_t j = ni & 63;
    const double n = s - sShifter;
    const double r = x - n;

    double p = sC5 * r;
    p = p + sC4;
    p = p * r;
    p = p + sC3;
    p = p * r;
    p = p + sC2;
    p = p * r;
    p = p + sC1;
    p = p * r;
    const double thi = sTab[2 * j];
    p = p + sTab[2 * j + 1];
    p = p * thi;
    p = p + thi;  /* p = Thi + Thi*(2^r - 1 + Tlo), in [1, 2) */

    if (x < sDenormX) {  /* x < -1022: denormal result, scale by 2^(M+100), then by 2^-100 */
        const uint32_t be = ((ni >> 6) + 1123) & 2047;
        p = bitsd((uint64_t)be << 52) * p;
        *pr = p * sTwoM100;
        return 4;
    }
    uint32_t be = ((ni >> 6) + 1023) & 2047;  /* biased exponent of 2^M */
    if (be <= 2046) {
        *pr = bitsd((uint64_t)be << 52) * p;
    }
    else {  /* M = 1024: scale by 2^1023, then by 2 */
        be = (be - 1) & 2047;
        p = bitsd((uint64_t)be << 52) * p;
        *pr = p * sTwo;
    }
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define RD (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_exp28_ha(__m512d x)
{
    /* r = x - floor(16 x)/16: imm 0x41 = 4 fraction bits, round down (also for the subtraction) */
    const __m512d r = _mm512_reduce_round_pd(x, 0x41, _MM_FROUND_NO_EXC);
    __m512d p = _mm512_fmadd_round_pd(V1(vC6), r, V1(vC5), RN);
    p = _mm512_fmadd_round_pd(r, p, V1(vC4), RN);
    /* t = x + 1.5*2^48 rounded down: its low 4 mantissa bits are floor(16 x) mod 16 */
    const __m512d t = _mm512_add_round_pd(x, V1(vShifter), RD);
    const __m512d ax = _mm512_and_pd(x, _mm512_castsi512_pd(_mm512_set1_epi64((long long)vAbsMaskBits)));
    p = _mm512_fmadd_round_pd(r, p, V1(vC3), RN);
    const __mmask8 callout = _mm512_cmp_round_pd_mask(ax, V1(vThreshold), _CMP_GE_OQ, _MM_FROUND_NO_EXC);
    p = _mm512_fmadd_round_pd(r, p, V1(vC2), RN);
    p = _mm512_fmadd_round_pd(r, p, V1(vC1), RN);

    const __m512i idx = _mm512_and_si512(_mm512_castpd_si512(t), _mm512_set1_epi64((long long)vIndexMaskBits));
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));
    const __m512d tl = _mm512_permutex2var_pd(_mm512_loadu_pd(vTl), idx, _mm512_loadu_pd(vTl + 8));

    const __m512d thr = _mm512_mul_round_pd(th, r, RN);         /* Th*r */
    p = _mm512_fmadd_round_pd(thr, p, tl, RN);                   /* Tl + Th*r*P(r) */
    const __m512d y = _mm512_add_round_pd(p, th, RN);            /* Th + Tl + Th*r*P(r) */
    __m512d res = _mm512_scalef_round_pd(y, x, RN);              /* * 2^floor(x) */

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dexp2_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_exp2_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_exp28_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_exp28_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
