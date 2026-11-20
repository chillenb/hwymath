/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision cosh, high-accuracy variant
 * (svml_z0_cosh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with a = |x|:
 *   m  = a*log2(e) + 1.5*2^48       (round to nearest: N = m - 1.5*2^48 is a*log2(e) rounded to a
 *                                    multiple of 1/16)
 *   N  = M + j/16, j = low 4 bits of m
 *   r  = a - N*ln2_hi - N*ln2_lo    (|r| <= ln2/32)
 *   Tp = 2^(j/16)/2  * 2^M = TpH + TpL     (16-entry tables, two-source permutes; the scaling by 2^M
 *   Tm = 2^(-j/16)/2 * 2^-M ~ TmH           is an integer addition to the exponent field)
 *   cosh(a) = (Tp + Tm) cosh(r) + (Tp - Tm) sinh(r)
 *           = TpH + [TmH + ((Tp - Tm)(r + s3 r^3 + s5 r^5 + s7 r^7)
 *                           + ((Tp + Tm)(c2 r^2 + c4 r^4 + c6 r^6) + TpL))]
 * Lanes with |x| >= 707.70 (high word above 0x40861d99; includes inf and NaN) are recomputed by a
 * scalar routine (exp-based for large |x|, signals overflow above 710.476).
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_cosh8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../cosh_d.h"

#if defined(__FAST_MATH__)
#error "svml_cosh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_cosh_d_ha_constants.inc"

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

/* 2^(e - 1023) for a biased exponent field e (0..2047), as the assembly builds it */
static inline double
pow2_field(uint32_t e)
{
    return bitsd((uint64_t)(e & 2047) << 52);
}

/*
 * __svml_dcosh_ha_cout_rare_internal: scalar cosh with a 64-entry table of 2^(j/64). Returns 0, or 3
 * (overflow). Statement order mirrors the assembly; FP flags are raised as by the assembly.
 *   inf, NaN:                 x*x
 *   |x| < 2^-53:              1 + |x|
 *   2^-53 <= |x| < 21.49:     (e^|x| + e^-|x|)/2 = Tp cosh(r) + ... with Tp = 2^(M-1) 2^(j/64),
 *                             Tm = 2^(-M-1) 2^(-j/64)
 *   21.49 <= |x| < 710.476:   e^|x|/2
 *   |x| >= 710.476:           +inf (overflow, return code 3)
 */
int
svml_dcosh_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t ex = (uint32_t)(xb >> 52) & 0x7ff;

    if (ex == 0x7ff) {  /* inf or NaN */
        *pr = x * x;
        return 0;
    }
    const double ax = bitsd(xb & 0x7fffffffffffffffULL);
    if (ex <= 969) {  /* |x| < 2^-53 */
        *pr = sOne + ax;
        return 0;
    }
    if (!(sOverflowX > ax)) {  /* |x| >= 710.476: +inf, raises overflow */
        volatile double h = sHuge;  /* computed at run time to raise the overflow flag */
        *pr = h * h;
        return 3;
    }

    if (!(ax < sLargeX)) {  /* 21.49 <= |x| < 710.476: e^|x| / 2 */
        double t = s64OverLn2 * ax;
        t = t + sShifter;
        const uint32_t ni = (uint32_t)dbits(t);
        const uint32_t j = ni & 63;
        const double n = t - sShifter;
        const double nh = n * sLn2o64Hi;
        const double nl = n * sLn2o64Lo;
        double r = ax - nh;
        r = r - nl;
        double p = sP6 * r;
        p = p + sP5;
        p = p * r;
        p = p + sP4;
        p = p * r;
        p = p + sP3;
        p = p * r;
        p = p + sP2;
        p = p * r;
        p = p * r;
        p = p + r;
        p = p * sTab[2 * j];
        p = p + sTab[2 * j + 1];
        p = p + sTab[2 * j];
        uint32_t be = ((ni >> 6) + 1022) & 2047;  /* biased exponent of 2^(M-1) */
        if (be <= 2046) {
            *pr = pow2_field(be) * p;
        }
        else {
            be = (be - 1) & 2047;
            p = p * pow2_field(be);
            *pr = p * sTwo;
        }
        return 0;
    }

    /* 2^-53 <= |x| < 21.49: |x| = (M*64 + j) ln2/64 + r,
     * cosh|x| = (Tp + Tm) cosh(r) + (Tp - Tm) sinh(r) */
    double t = s64OverLn2 * ax;
    t = t + sShifter;
    const uint32_t ni = (uint32_t)dbits(t);
    const uint32_t m = ni >> 6;
    const uint32_t j = ni & 63;
    const double n = t - sShifter;
    const double scale_p = pow2_field(m + 1022);                       /* 2^(M-1) */
    const double nh = n * sLn2o64Hi;
    const double scale_m = pow2_field((uint32_t)(-(m + 1023)) - 4);    /* 2^(-M-2) */
    const double nl = n * sLn2o64Lo;
    double r = ax - nh;
    r = r - nl;
    const double tph = sTab[2 * j] * scale_p;
    const double r2 = r * r;
    const double tmh = sTab[128 - 2 * j] * scale_m;  /* 2^(1-j/64) 2^(-M-2) */
    double odd = sP5 * r2;
    const double s = tph + tmh;
    const double d = tph - tmh;
    odd = odd + sP3;
    const double tpl = sTab[2 * j + 1] * scale_p;
    double sl = tph - s;
    const double tml = sTab[129 - 2 * j] * scale_m;
    odd = odd * r2;
    sl = sl + tmh;
    odd = odd * r;  /* sinh(r) - r */
    double even = sP6 * r2;
    double lo = tpl - tml;
    odd = odd * d;
    lo = lo * r;
    even = even + sP4;
    const double rd = r * d;
    even = even * r2;
    lo = tpl + lo;
    lo = lo + tml;
    even = even + sP2;
    lo = lo + odd;
    even = even * r2;  /* cosh(r) - 1 */
    lo = lo + rd;
    even = even * s;
    lo = lo + even;
    lo = lo + sl;
    *pr = s + lo;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define VB(b) _mm512_castsi512_pd(_mm512_set1_epi64((long long)(b)))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_cosh8_ha(__m512d x)
{
    const __m512d S = V1(vShifter);
    const __m512d ax = _mm512_andnot_pd(VB(vSignMaskBits), x);
    const __m512d m = _mm512_fmadd_round_pd(ax, V1(vLog2e), S, RN);  /* a*log2(e) + 1.5*2^48 */
    const __m512i mb = _mm512_castpd_si512(m);  /* vpermt2pd uses its low 4 bits, j */

    /* callout lanes: high word of |x| > 0x40861d99 (signed 32-bit compare; inf and NaN included) */
    const __m256i hi = _mm512_cvtepi64_epi32(_mm512_srli_epi64(_mm512_castpd_si512(ax), 32));
    const __m512d n = _mm512_sub_round_pd(m, S, RN);
    const __m512i tpl = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTpL), mb, _mm512_loadu_pd(vTpL + 8)));
    const __m512i tmh = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTmH), mb, _mm512_loadu_pd(vTmH + 8)));
    __m512d r = _mm512_fnmadd_round_pd(n, V1(vLn2Hi), ax, RN);  /* a - N*ln2_hi */
    r = _mm512_fnmadd_round_pd(n, V1(vLn2Lo), r, RN);          /* - N*ln2_lo */
    const __m512i m48 = _mm512_slli_epi64(mb, 48);
    const __m256i big = _mm256_cmpgt_epi32(hi, _mm256_set1_epi32((int)vThresholdHi));
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const int callout = _mm256_movemask_ps(_mm256_castsi256_ps(big));
    const __m512i tph = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTpH), mb, _mm512_loadu_pd(vTpH + 8)));

    /* scaling by 2^+-M: M << 52 added to / subtracted from the bit patterns */
    const __m512i e = _mm512_and_epi64(m48, _mm512_set1_epi64((long long)vExpMaskBits));
    const __m512d tpl2 = _mm512_castsi512_pd(_mm512_add_epi64(e, tpl));
    const __m512d tp = _mm512_castsi512_pd(_mm512_add_epi64(e, tph));
    const __m512d tm = _mm512_castsi512_pd(_mm512_sub_epi64(tmh, e));

    const __m512d sum = _mm512_add_round_pd(tp, tm, RN);
    const __m512d dif = _mm512_sub_round_pd(tp, tm, RN);
    __m512d po = _mm512_fmadd_round_pd(r2, V1(vS7), V1(vS5), RN);
    po = _mm512_fmadd_round_pd(r2, po, V1(vS3), RN);
    __m512d pe = _mm512_fmadd_round_pd(r2, V1(vC6), V1(vC4), RN);
    po = _mm512_mul_round_pd(po, r2, RN);                    /* r^2 (s3 + s5 r^2 + s7 r^4) */
    pe = _mm512_fmadd_round_pd(r2, pe, V1(vC2), RN);
    po = _mm512_fmadd_round_pd(r, po, r, RN);                /* r + r^3 (...) */
    pe = _mm512_mul_round_pd(pe, r2, RN);                    /* r^2 (c2 + c4 r^2 + c6 r^4) */
    pe = _mm512_fmadd_round_pd(sum, pe, tpl2, RN);           /* (Tp + Tm)(cosh(r) - 1) + TpL */
    __m512d res = _mm512_fmadd_round_pd(dif, po, pe, RN);    /* + (Tp - Tm) sinh(r) */
    res = _mm512_add_round_pd(res, tm, RN);
    res = _mm512_add_round_pd(res, tp, RN);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1 << i)) {
                svml_dcosh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_cosh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_cosh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 k = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, k, svml_cosh8_ha(_mm512_maskz_loadu_pd(k, x + i)));
    }
}

#endif /* __x86_64__ */
