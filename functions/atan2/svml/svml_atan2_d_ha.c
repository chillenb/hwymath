/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision atan2, high-accuracy variant
 * (svml_z0_atan2_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes; ax = |x|, ay = |y|): the interval of |y/x| is found without a
 * division, by comparing ay with ax*7/16, ax*11/16, ax*19/16 and ax*39/16:
 *   |y/x| <  7/16:  a = ay,              b = ax,              atan(0)   = 0
 *   |y/x| < 11/16:  a = ay - 0.5 ax,     b = ax + 0.5 ay,     atan(1/2) = Hi + Lo
 *   |y/x| < 19/16:  a = ay - ax,         b = ax + ay,         atan(1)   = pi/4
 *   |y/x| < 39/16:  a = ay - 1.5 ax,     b = ax + 1.5 ay,     atan(3/2)
 *   otherwise:      a = -ax,             b = ay,              atan(inf) = pi/2
 * (one FMA each), so that atan(ay/ax) = atan(c) + atan(a/b) with |a/b| <= 7/16. The quotient
 *   r0 = rcp14(b) (vrcp14pd, |rel err| < 2^-14), refined twice: r1 = r0 + r0(1 - b r0), r2 = r1 + r1(1 - b r1)
 *   R = RN(r2 a), Rlo = r2 (a - R b)                      (R + Rlo = a/b in double-double)
 *   atan(a/b) = R + (R*R2*P(R2) + Rlo),  P of degree 11 in R2 = R^2 (even/odd halves in R^4)
 * The low parts (Rlo, atan(c)_lo, and for x < 0 the low part of pi, with the sign of x) are added
 * first, then R, then atan(c)_hi; the result gets the sign of x, pi is added for x < 0 (atan2 =
 * pi - atan(ay/ax)), and the sign of y is or-ed in.
 * Lanes where x or y has a high word outside [0x00300000, 0x7e000000) (|v| < 2^-1020, |v| >= 2^993, inf,
 * NaN) are special: those with a zero argument (and no NaN) are computed in the vector code
 * (+-0, +-pi/2, +-pi), the others by the scalar routine.
 *
 * The scalar routine (__svml_datan2_ha_cout_rare_internal: a complete scalar atan2 with its own
 * table-driven algorithm, scaling for extreme exponents and the special cases) is translated
 * mechanically, instruction by instruction, by functions/asin/tools/asm2c.py into
 * svml_datan2_ha_cout_rare.inc.
 *
 * vrcp14pd results are implementation-defined (Intel's are a function of the top 16 mantissa bits,
 * with at most 16 significant bits) and the result depends on them (through the rounding of r2 and R).
 * Intel publishes a bit-exact reference
 * (RECIP14.c), which the Highway version emulates (functions/asin/hwy/approx14-inl.h).
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only the vector functions need AVX-512F at run time; the scalar routines are plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../atan2_d.h"
#include "asm_model.h"

#if defined(__FAST_MATH__)
#error "svml_atan2_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routines use separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_atan2_d_ha_constants.inc"

/* __svml_datan2_ha_cout_rare_internal */
#include "svml_datan2_ha_cout_rare.inc"

/* The vector code's zero-argument branch (.LBL_1_12) for one lane with x or y = +-0 and no NaN:
 * PIO2 = 0 if |y| < |x| else pi/2, set to 0 if max(|x|, |y|) = 0; (sign(x) | PIO2) + pi if x < 0
 * (sign bit), and the sign of y or-ed in. */
double
svml_atan2_ha_zero_args(double y, double x)
{
    const double ax = bitsd(dbits(x) & vAbsMaskBits), ay = bitsd(dbits(y) & vAbsMaskBits);
    const int k1 = ay < ax;
    const double den = k1 ? ax : ay;
    const uint64_t pio2 = k1 ? 0 : dbits(vZPio2);
    const uint64_t v = den != 0.0 ? pio2 : dbits(vZero);
    double r = bitsd((dbits(x) & vSignMaskBits) | v);
    if (dbits(x) >> 63) r = r + vZPi;
    return bitsd(dbits(r) | (dbits(y) & vSignMaskBits));
}

#if defined(__x86_64__)  /* the vector routines; the scalar routines above are portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define AVX512 __attribute__((target("avx512f")))
#define PD(v) _mm512_castsi512_pd(v)
#define SI(v) _mm512_castpd_si512(v)
/* compiler barrier without instructions: keeps the instruction mix comparable with the assembly */
#define OPAQUE(v) __asm__("" : "+v"(v))

AVX512 static inline __attribute__((always_inline)) __m512d
atan2_8(__m512d y, __m512d x)
{
    const __m512i absmask = _mm512_set1_epi64((long long)vAbsMaskBits);
    __m512d ax = PD(_mm512_and_si512(SI(x), absmask));
    __m512d ay = PD(_mm512_and_si512(SI(y), absmask));
    OPAQUE(ax);
    OPAQUE(ay);
    /* 0.0 read from memory, as in the assembly (a vxorpd zeroing idiom would count as an extra xor) */
    const __m512d zero = _mm512_set1_pd(asm_model_zero);
    const __mmask8 xneg = _mm512_cmp_round_pd_mask(x, zero, _CMP_LT_OQ, SAE);

    /* special lanes: high word of |x| or |y| outside [0x00300000, 0x7e000000) */
    const __m256i hbase = _mm256_set1_epi32((int)vHiBase), hlim = _mm256_set1_epi32((int)vHiLimit);
    const __m256i hx = _mm256_sub_epi32(_mm512_cvtepi64_epi32(_mm512_srli_epi64(SI(ax), 32)), hbase);
    const __m256i hy = _mm256_sub_epi32(_mm512_cvtepi64_epi32(_mm512_srli_epi64(SI(ay), 32)), hbase);
    __m512d sgnx = PD(_mm512_xor_si512(SI(x), SI(ax)));
    __m512d sgny = PD(_mm512_xor_si512(SI(y), SI(ay)));
    OPAQUE(sgnx);
    OPAQUE(sgny);
    const __m256i spx = _mm256_or_si256(_mm256_cmpgt_epi32(hx, hlim), _mm256_cmpeq_epi32(hx, hlim));
    const __m256i spy = _mm256_or_si256(_mm256_cmpgt_epi32(hy, hlim), _mm256_cmpeq_epi32(hy, hlim));

    /* interval of |y/x| */
    const __m512d t7 = _mm512_mul_round_pd(ax, V1(vT7_16), RN);
    const __m512d t11 = _mm512_mul_round_pd(ax, V1(vT11_16), RN);
    const __m512d t19 = _mm512_mul_round_pd(ax, V1(vT19_16), RN);
    const __mmask8 k1 = _mm512_cmp_round_pd_mask(t11, ay, _CMP_LT_OQ, SAE);
    const __mmask8 k5 = _mm512_cmp_round_pd_mask(t7, ay, _CMP_LT_OQ, SAE);
    const __mmask8 k2 = _mm512_cmp_round_pd_mask(t19, ay, _CMP_LT_OQ, SAE);
    __m256i special = _mm256_or_si256(spy, spx);
    __m512d hi = _mm512_mask_blend_pd(k1, V1(vAtan05Hi), V1(vPio4Hi));
    const __m512d t39 = _mm512_mul_round_pd(ax, V1(vT39_16), RN);
    const __mmask8 k3 = _mm512_cmp_round_pd_mask(t39, ay, _CMP_LT_OQ, SAE);
    __m512d d = _mm512_mask_blend_pd(k1, V1(vHalf), V1(vOne));
    __m512d lo = _mm512_mask_blend_pd(k1, V1(vAtan05Lo), V1(vPio4Lo));
    const __m512d d2 = _mm512_mask_blend_pd(k3, V1(vOneHalf), V1(vOne));
    const __m512d hi2 = _mm512_mask_blend_pd(k3, V1(vAtan15Hi), V1(vPio2Hi));
    d = _mm512_mask_blend_pd(k2, d, d2);
    const __m512d lo2 = _mm512_mask_blend_pd(k3, V1(vAtan15Lo), V1(vPio2Lo));
    hi = _mm512_mask_blend_pd(k2, hi, hi2);
    lo = _mm512_mask_blend_pd(k2, lo, lo2);
    /* b = c*ax + d*ay (c = 0 in the last interval) */
    __m512d b = PD(_mm512_mask_xor_epi64(SI(ax), k3, SI(ax), SI(ax)));
    b = _mm512_mask3_fmadd_round_pd(d, ay, b, k5, RN);

    /* 1/b ~ r2, by two refinements of r0 */
    __m512d r = _mm512_rcp14_pd(b);
    const __m512d e = _mm512_fnmadd_round_pd(b, r, V1(vOne), RN);
    r = _mm512_fmadd_round_pd(e, r, r, RN);
    const __m512d e1 = _mm512_fnmadd_round_pd(b, r, V1(vOne), RN);
    /* a = c*ay - d*ax */
    __m512d a = PD(_mm512_mask_xor_epi64(SI(ay), k3, SI(ay), SI(ay)));
    r = _mm512_fmadd_round_pd(e1, r, r, RN);
    a = _mm512_mask3_fnmadd_round_pd(d, ax, a, k5, RN);

    /* R + Rlo = a/b */
    const __m512d rr = _mm512_mul_round_pd(r, a, RN);
    const __m512d r2 = _mm512_mul_round_pd(rr, rr, RN);
    const __m512d t = _mm512_fnmadd_round_pd(rr, b, a, RN);  /* a - R*b */
    const __m512d r4 = _mm512_mul_round_pd(r2, r2, RN);
    __m512d rlo = _mm512_mul_round_pd(t, r, RN);
    rlo = _mm512_mask_add_round_pd(rlo, k5, rlo, lo, RN);

    /* polynomial */
    __m512d pa = _mm512_fmadd_round_pd(V1(vC12), r4, V1(vC10), RN);
    __m512d pb = _mm512_fmadd_round_pd(V1(vC11), r4, V1(vC9), RN);
    pa = _mm512_fmadd_round_pd(r4, pa, V1(vC8), RN);
    pb = _mm512_fmadd_round_pd(r4, pb, V1(vC7), RN);
    pa = _mm512_fmadd_round_pd(r4, pa, V1(vC6), RN);
    pb = _mm512_fmadd_round_pd(r4, pb, V1(vC5), RN);
    pa = _mm512_fmadd_round_pd(r4, pa, V1(vC4), RN);
    pb = _mm512_fmadd_round_pd(r4, pb, V1(vC3), RN);
    pa = _mm512_fmadd_round_pd(r4, pa, V1(vC2), RN);
    pb = _mm512_fmadd_round_pd(r4, pb, V1(vC1), RN);
    __m512d p = _mm512_fmadd_round_pd(r2, pa, pb, RN);
    p = _mm512_mul_round_pd(p, r2, RN);

    /* low parts, then R, then atan(c)_hi; sign of x, + pi for x < 0, sign of y */
    const __m512d pilo = PD(_mm512_xor_si512(SI(sgnx), _mm512_set1_epi64((long long)dbits(vPiLo))));
    rlo = _mm512_mask_add_round_pd(rlo, xneg, rlo, pilo, RN);
    p = _mm512_fmadd_round_pd(rr, p, rlo, RN);
    __m512d res = _mm512_add_round_pd(rr, p, RN);
    res = _mm512_mask_add_round_pd(res, k5, res, hi, RN);
    res = PD(_mm512_xor_si512(SI(res), SI(sgnx)));
    res = _mm512_mask_add_round_pd(res, xneg, res, V1(vPiHi), RN);
    res = PD(_mm512_or_si512(SI(res), SI(sgny)));
    OPAQUE(res);

    unsigned callout = 0;
    if (_mm256_movemask_ps(_mm256_castsi256_ps(special))) {
        /* zero-argument branch (.LBL_1_12) */
        const __m512i ones = _mm512_set1_epi64(-1);
        const __mmask8 xnan = _mm512_cmp_round_pd_mask(x, x, _CMP_UNORD_Q, SAE);
        const __mmask8 ynan = _mm512_cmp_round_pd_mask(y, y, _CMP_UNORD_Q, SAE);
        const __mmask8 ylt = _mm512_cmp_round_pd_mask(ay, ax, _CMP_LT_OQ, SAE);  /* PIO2 = 0 if |y| < |x| */
        const __m512d den = _mm512_mask_blend_pd(ylt, ay, ax);
        const __m512i signmask = _mm512_set1_epi64((long long)vSignMaskBits);
        const __m512i sx = _mm512_and_si512(SI(x), signmask), sy = _mm512_and_si512(SI(y), signmask);
        __m512i pio2 = _mm512_set1_epi64((long long)dbits(vZPio2));
        pio2 = _mm512_mask_xor_epi64(pio2, ylt, pio2, pio2);
        const __m512i xok = _mm512_mask_andnot_epi64(ones, xnan, SI(x), SI(x));
        const __m512i yok = _mm512_mask_andnot_epi64(ones, ynan, SI(y), SI(y));
        const __m512i notnan = _mm512_srli_epi64(_mm512_and_si512(xok, yok), 32);
        const __mmask8 xnz = _mm512_cmp_round_pd_mask(ax, V1(vZero), _CMP_NEQ_UQ, SAE);
        const __mmask8 ynz = _mm512_cmp_round_pd_mask(ay, V1(vZero), _CMP_NEQ_UQ, SAE);
        const __mmask8 dnz = _mm512_cmp_round_pd_mask(den, V1(vZero), _CMP_NEQ_UQ, SAE);
        const __mmask8 xsign = _mm512_cmpgt_epi64_mask(SI(V1(vZero)), SI(x));
        const __m256i notnan32 = _mm512_cvtepi64_epi32(notnan);
        const __m512i xz = _mm512_mask_andnot_epi64(ones, xnz, SI(ax), SI(ax));
        const __m512i yz = _mm512_mask_andnot_epi64(ones, ynz, SI(ay), SI(ay));
        const __m512i anyz = _mm512_or_si512(yz, xz);
        const __m256i zero_ok = _mm256_and_si256(notnan32, _mm512_cvtepi64_epi32(_mm512_srli_epi64(anyz, 32)));
        special = _mm256_andnot_si256(zero_ok, special);  /* zero (and no NaN) arguments: no callout */
        const __m512i zero_ok64 = _mm512_cvtepu32_epi64(zero_ok);
        callout = (unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(special));
        const __m512i denz = _mm512_mask_andnot_epi64(ones, dnz, SI(den), SI(den));
        const __m512i v = _mm512_or_si512(_mm512_and_si512(denz, SI(V1(vZero))), _mm512_andnot_si512(denz, pio2));
        __m512d z = PD(_mm512_or_si512(sx, v));
        z = _mm512_mask_add_round_pd(z, xsign, z, V1(vZPi), RN);
        const __m512i zres = _mm512_or_si512(SI(z), sy);
        const __m512i sel = _mm512_or_epi32(zero_ok64, _mm512_slli_epi64(zero_ok64, 32));
        res = PD(_mm512_or_si512(_mm512_and_si512(sel, zres), _mm512_andnot_si512(sel, SI(res))));
    }
    if (callout) {
        double yin[8], xin[8], out[8];
        _mm512_storeu_pd(yin, y);
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_datan2_ha_cout_rare(&yin[i], &xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

AVX512 __m512d
svml_atan28_ha(__m512d y, __m512d x)
{
    return atan2_8(y, x);
}

#define ARRAY(NAME, KERNEL)                                                                         \
    AVX512 void NAME(const double *y, const double *x, double *r, size_t len)                      \
    {                                                                                               \
        size_t i = 0;                                                                               \
        for (; i + 8 <= len; i += 8) {                                                              \
            _mm512_storeu_pd(r + i, KERNEL(_mm512_loadu_pd(y + i), _mm512_loadu_pd(x + i)));        \
        }                                                                                           \
        if (i < len) {                                                                              \
            const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);                                   \
            _mm512_mask_storeu_pd(r + i, m, KERNEL(_mm512_maskz_loadu_pd(m, y + i), _mm512_maskz_loadu_pd(m, x + i))); \
        }                                                                                           \
    }
ARRAY(svml_atan2_ha_array, svml_atan28_ha)

#endif /* __x86_64__ */
