/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision sinh, high-accuracy variant
 * (svml_z0_sinh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with s = sign(x), a = |x|:
 *   m  = a*log2(e) + 1.5*2^48       (round to nearest: N = m - 1.5*2^48 is a*log2(e) rounded to a
 *                                    multiple of 1/16; if N = 1/16 exactly, N = 0 is used instead)
 *   N  = M + j/16, j = low 4 bits of m
 *   r  = a - N*ln2_hi - N*ln2_lo    (|r| <= ln2/32)
 *   Tp = 2^(j/16)/2  * 2^M = TpH + TpL     (16-entry tables, two-source permutes; the scaling by 2^M
 *   Tm = 2^(-j/16)/2 * 2^-M = TmH + TmL     is an integer addition to the exponent field; TmL is set
 *                                            to 0 where that would underflow)
 *   sinh(a) = (Tp - Tm) cosh(r) + (Tp + Tm) sinh(r)
 *           = D + [D*(c2 r^2 + ... + c8 r^8) + lo] + S*(r + s3 r^3 + ... + s7 r^7)
 *   with D = TpH - TmH, S = TpH + TmH, lo = (TpL - (rounding error of D)) - TmL;
 *   result = s * sinh(a).
 * Lanes with |x| >= 669.78 (high word above 0x4084ee33; includes inf and NaN) are recomputed by a
 * scalar routine (exp-based for large |x|, signals overflow above 710.476).
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_sinh8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../sinh_d.h"

#if defined(__FAST_MATH__)
#error "svml_sinh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_sinh_d_ha_constants.inc"

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

/* |v| with the sign of x (the assembly overwrites the sign bit of the result's top byte) */
static inline double
with_sign_of(double v, uint64_t xb)
{
    return bitsd((dbits(v) & 0x7fffffffffffffffULL) | (xb & 0x8000000000000000ULL));
}

/* 2^(e - 1023) for a biased exponent field e (0..2047), as the assembly builds it */
static inline double
pow2_field(uint32_t e)
{
    return bitsd((uint64_t)(e & 2047) << 52);
}

/*
 * __svml_dsinh_ha_cout_rare_internal: scalar sinh with a 64-entry table of 2^(j/64). Returns 0, or 3
 * (overflow). Statement order mirrors the assembly; FP flags are raised as by the assembly.
 *   |x| < 2^-53 (incl. denormals and 0): x (raising inexact/underflow like the assembly)
 *   |x| < 0.0272:                         x + x^3 (s3 + s5 x^2 + ...), Horner
 *   0.0272 <= |x| < 21.49:                (e^|x| - e^-|x|)/2 in double-double arithmetic
 *   21.49 <= |x| < 710.476:               e^|x|/2
 *   |x| >= 710.476:                       +-inf (overflow, return code 3)
 */
int
svml_dsinh_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t ex = (uint32_t)(xb >> 52) & 0x7ff;

    if (ex == 0x7ff) {  /* inf or NaN */
        *pr = x + x;
        return 0;
    }
    if (ex == 0) {  /* zero or denormal: x*tiny + x */
        volatile double t = x * sTiny;
        *pr = t + x;
        return 0;
    }
    const double ax = bitsd(xb & 0x7fffffffffffffffULL);
    if (ex <= 969) {  /* |x| < 2^-53 */
        volatile double one = sTab[0];
        const double c = one + sTiny;  /* 1.0, raising inexact */
        *pr = c * x;
        return 0;
    }
    if (!(sOverflowX > ax)) {  /* |x| >= 710.476: +-inf, raises overflow */
        *pr = sHuge * x;
        return 3;
    }

    if (!(sLargeX > ax)) {  /* 21.49 <= |x| < 710.476: e^|x| / 2 */
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
        double res;
        if (be <= 2046) {
            res = p * pow2_field(be);
        }
        else {
            be = (be - 1) & 2047;
            res = p * pow2_field(be);
            res = res * sTwo;
        }
        *pr = with_sign_of(res, xb);
        return 0;
    }

    if (ax < sSmallX) {  /* 2^-53 <= |x| < 0.0272: odd Taylor polynomial */
        double x2 = ax * ax;
        double p = sS9 * x2;
        p = p + sS7;
        p = p * x2;
        p = p + sS5;
        p = p * x2;
        p = p + sS3;
        x2 = x2 * p;
        x2 = x2 * ax;
        const double res = ax + x2;
        *pr = with_sign_of(res, xb);
        return 0;
    }

    /* 0.0272 <= |x| < 21.49: |x| = (M*64 + j) ln2/64 + r,
     * sinh|x| = (Tp - Tm) cosh(r) + (Tp + Tm) sinh(r), Tp = 2^(M-1) 2^(j/64), Tm = 2^(-M-1) 2^(-j/64) */
    double t = s64OverLn2 * ax;
    t = t + sShifter;
    const uint32_t ni = (uint32_t)dbits(t);
    const uint32_t m = ni >> 6;
    const uint32_t j = ni & 63;
    const double n = t - sShifter;
    const double nh = sLn2o64Hi * n;
    const double mnl = -n * sLn2o64Lo;
    const double a = ax - nh;
    /* r = rh + rl (double-double) */
    const double rh = mnl + a;
    const double d1 = a - rh;
    const double s48 = rh + d1;
    const double rl0 = mnl + d1;
    const double scale_p = pow2_field(m + 1022);               /* 2^(M-1) */
    const double tph = sTab[2 * j] * scale_p;
    const double tpl = sTab[2 * j + 1] * scale_p;
    const double scale_m = pow2_field((uint32_t)(-(m + 1023)) - 4);  /* 2^(-M-2) */
    const double a2 = a - s48;
    const double rl = rl0 + a2;
    /* split rh = rhh + rhl (Dekker) */
    const double sph = sSplitter * rh;
    const double r2 = rh * rh;
    const double sph2 = sph - rh;
    const double rhh = sph - sph2;
    double odd = sP5 * r2;   /* sinh(r) - r ~ rh r2 (1/6 + r2/120) */
    odd = odd + sP3;
    double even = sP6 * r2;  /* cosh(r) - 1 ~ r2 (1/2 + r2/24 + r2^2/720) */
    odd = odd * r2;
    even = even + sP4;
    odd = odd * rh;
    even = even * r2;
    even = even + sP2;
    const double rhl = rh - rhh;
    even = even * r2;
    /* D = Tp - Tm, S = Tp + Tm as double-doubles */
    const double tmh = sTab[128 - 2 * j] * scale_m;  /* 2^(1-j/64) 2^(-M-2) */
    const double d = tph - tmh;
    const double tml = sTab[129 - 2 * j] * scale_m;
    double e1 = tph - d;
    e1 = e1 - tmh;
    e1 = e1 + tpl;
    const double tl = tpl + tml;
    e1 = e1 - tml;
    const double dh = d + e1;
    const double s = tph + tmh;
    const double dd = d - dh;
    const double dl = e1 + dd;
    const double s2 = tph - s;
    double sl = tmh + s2;
    sl = sl + tl;
    const double sh = s + sl;
    const double s3 = s - sh;
    sl = sl + s3;
    const double sps = sSplitter * sh;
    const double rh_sl = rh * sl;
    const double sps2 = sps - sh;
    const double shh = sps - sps2;
    const double sh_odd = sh * odd;
    const double shl = sh - shh;
    const double sh_rl = sh * rl;
    const double rl_sl = rl * sl;
    const double odd_sl = odd * sl;
    const double dh_even = dh * even;
    double lo = even * dl;
    const double q = sh_odd + dh_even;
    lo = lo + dl;
    lo = lo + sh_rl;
    lo = lo + rl_sl;
    lo = lo + odd_sl;
    lo = lo + rh_sl;
    const double q2 = sh_odd - q;
    double ql = dh_even + q2;
    const double rhh_shh = rhh * shh;
    const double shh_rhl = shh * rhl;
    const double rhl_shl = rhl * shl;
    const double rhh_shl = rhh * shl;
    const double v2 = rhh_shh + q;
    ql = ql + shh_rhl;
    ql = ql + rhl_shl;
    const double v3 = rhh_shh - v2;
    ql = ql + rhh_shl;
    lo = lo + ql;
    const double v2l = q + v3;
    const double v4 = dh + v2;
    lo = lo + v2l;
    const double v5 = dh - v4;
    const double v4l = v2 + v5;
    lo = lo + v4l;
    const double res = v4 + lo;
    *pr = with_sign_of(res, xb);
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define VB(b) _mm512_castsi512_pd(_mm512_set1_epi64((long long)(b)))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
/* Optimization barrier: keeps the compiler from merging the bitwise operations into different ones
 * (andn/ternlog), so that the instruction mix stays that of the assembly (make sinh-insns). No effect
 * on results. */
#define OPAQUE(v) __asm__("" : "+v"(v))

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_sinh8_ha(__m512d x)
{
    const __m512d signmask = VB(vSignMaskBits);
    const __m512d S = V1(vShifter);
    __m512d sign = _mm512_and_pd(signmask, x);
    OPAQUE(sign);
    const __m512d ax = _mm512_xor_pd(sign, x);
    __m512d m = _mm512_fmadd_round_pd(ax, V1(vLog2e), S, RN);  /* a*log2(e) + 1.5*2^48 */

    /* callout lanes: high word of |x| > 0x4084ee33 (signed 32-bit compare; inf and NaN included) */
    const __m256i hi = _mm512_cvtepi64_epi32(_mm512_srli_epi64(_mm512_castpd_si512(ax), 32));
    const __mmask8 k1 = _mm512_cmp_round_pd_mask(VB(vShifterP1Bits), m, _CMP_EQ_OQ, _MM_FROUND_NO_EXC);
    const __m256i big = _mm256_cmpgt_epi32(hi, _mm256_set1_epi32((int)vThresholdHi));
    m = _mm512_mask_blend_pd(k1, m, S);  /* N = 1/16 -> N = 0 */
    const int callout = _mm256_movemask_ps(_mm256_castsi256_ps(big));

    const __m512i mb = _mm512_castpd_si512(m);
    __m512i j = _mm512_and_epi64(mb, _mm512_set1_epi64((long long)vIndexMaskBits));
    OPAQUE(j);
    const __m512i mhi = _mm512_xor_epi64(j, mb);  /* m without the index bits */
    __m512i tph = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTpH), j, _mm512_loadu_pd(vTpH + 8)));
    __m512i tmh = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTmH), j, _mm512_loadu_pd(vTmH + 8)));
    __m512i tpl = _mm512_castpd_si512(_mm512_permutex2var_pd(_mm512_loadu_pd(vTpL), j, _mm512_loadu_pd(vTpL + 8)));
    __m512d tml = _mm512_permutex2var_pd(_mm512_loadu_pd(vTmL), j, _mm512_loadu_pd(vTmL + 8));

    const __m512d n = _mm512_sub_round_pd(m, S, RN);
    const __m512d atml = _mm512_andnot_pd(signmask, tml);  /* |TmL| */
    __m512d r = _mm512_fnmadd_round_pd(n, V1(vLn2Hi), ax, RN);  /* a - N*ln2_hi */
    r = _mm512_fnmadd_round_pd(n, V1(vLn2Lo), r, RN);          /* - N*ln2_lo */

    /* scaling by 2^+-M: M << 52 added to / subtracted from the bit patterns */
    const __m512i e = _mm512_slli_epi64(mhi, 48);
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __mmask8 k2 = _mm512_cmp_round_pd_mask(atml, _mm512_castsi512_pd(e), _CMP_LT_OQ, _MM_FROUND_NO_EXC);
    const __m512d tp = _mm512_castsi512_pd(_mm512_add_epi64(e, tph));
    const __m512d tm = _mm512_castsi512_pd(_mm512_sub_epi64(tmh, e));
    const __m512d tpl2 = _mm512_castsi512_pd(_mm512_add_epi64(e, tpl));
    tml = _mm512_castsi512_pd(_mm512_sub_epi64(_mm512_castpd_si512(tml), e));
    const __m512d dif = _mm512_sub_round_pd(tp, tm, RN);
    const __m512d sum = _mm512_add_round_pd(tp, tm, RN);
    tml = _mm512_mask_xor_pd(tml, k2, tml, tml);  /* TmL*2^-M would underflow: 0 */
    __m512d c = _mm512_sub_round_pd(dif, tp, RN);
    c = _mm512_add_round_pd(tm, c, RN);  /* -(rounding error of dif) */

    /* poly = dif*(c2 r^2 + ... + c8 r^8) + lo + sum*(r + s3 r^3 + s5 r^5 + s7 r^7) */
    __m512d po = _mm512_fmadd_round_pd(r2, V1(vS7), V1(vS5), RN);
    __m512d lo = _mm512_sub_round_pd(tpl2, c, RN);
    lo = _mm512_sub_round_pd(lo, tml, RN);
    po = _mm512_fmadd_round_pd(r2, po, V1(vS3), RN);
    __m512d pe = _mm512_fmadd_round_pd(r2, V1(vC8), V1(vC6), RN);
    po = _mm512_mul_round_pd(po, r2, RN);
    pe = _mm512_fmadd_round_pd(r2, pe, V1(vC4), RN);
    po = _mm512_fmadd_round_pd(r, po, r, RN);
    pe = _mm512_fmadd_round_pd(r2, pe, V1(vC2), RN);
    pe = _mm512_mul_round_pd(r2, pe, RN);
    pe = _mm512_fmadd_round_pd(dif, pe, lo, RN);
    __m512d res = _mm512_fmadd_round_pd(sum, po, pe, RN);
    res = _mm512_add_round_pd(res, dif, RN);
    res = _mm512_or_pd(sign, res);

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1 << i)) {
                svml_dsinh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_sinh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_sinh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 k = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, k, svml_sinh8_ha(_mm512_maskz_loadu_pd(k, x + i)));
    }
}

#endif /* __x86_64__ */
