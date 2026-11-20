/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision sin, high-accuracy variant
 * (svml_z0_sin_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Main path (|x| <= 2^24), 8 lanes:
 *   y  = x*16/pi + 1.5*2^52 (round to nearest): N = y - 1.5*2^52 = round(x*16/pi), j = N mod 16
 *   R + Rlo = x - N*(PI1 + PI2 + PI3)          (pi/16 in three parts; Rlo collects the rounding errors)
 *   B = N*pi/16 (mod pi, j = N mod 16): Sh + Sl = sin B, Ch + Cl = cos B from 16-entry tables
 *     (vpermt2pd; Cl and Sl share one 64-bit table entry: Cl is the whole word, Sl its low half << 32)
 *   sin(B + R) = Sh + R*Ch [split exactly] + R^2*(Sh*Cp(R^2) + R*Ch*Sp(R^2))
 *                + Sl + Rlo*(Ch - Sh*R) + R*Cl + (error of R*Ch)
 *     with cos R - 1 = R^2*Cp, sin R - R = R^3*Sp (degree-3 polynomials in R^2);
 *   bit 4 of N (N mod 32 >= 16: B in [pi, 2pi)) flips the sign. For x = -0, N's bit 4 is flipped so that
 *   sin(-0) = -0.
 * Large path (|x| > 2^24, +-inf, NaN; all lanes recomputed if any lane needs it): Payne-Hanek reduction
 *   of |x| modulo pi/256 with a 2048x192-bit table of 2/pi (vpmullq, gathers), evaluation with a
 *   512-entry table of sin/cos(M*pi/256); see svml_trig_reduction.h. The sign of x is applied at the end.
 * Lanes with x = +-inf or NaN are then recomputed by the scalar callout (sin(+-inf) = NaN, invalid).
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_sin8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../sin_d.h"

#if defined(__FAST_MATH__)
#error "svml_sin_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_sin_d_ha_constants.inc"

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
 * __svml_dsin_ha_cout_rare_internal: called for the lanes with x = +-inf or NaN. sin(+-inf) = 0*x =
 * NaN with the invalid flag (return code 1), sin(NaN) = x*x (quiet NaN, invalid for signaling NaNs).
 * Other inputs: returns 0 without writing *pr. The x86 result of 0*inf is the default NaN
 * 0xfff8000000000000; it is stored explicitly so that every architecture gives the same bits.
 */
int
svml_dsin_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    if (((xb >> 52) & 0x7ff) != 0x7ff) {
        return 0;
    }
    if ((uint32_t)xb == 0 && (uint32_t)(xb >> 32 & 0x7fffffff) == 0x7ff00000) {  /* +-inf */
        volatile double zero = 0.0;  /* __dsin_ha_CoutTab[0]; multiplied at run time: raises invalid */
        volatile double r = x * zero;
        (void)r;
        *pr = bitsd(0xfff8000000000000ULL);
        return 1;
    }
    *pr = x * x;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#include "svml_trig_reduction.h"

#define V1(c) _mm512_set1_pd(c)
#define U1(c) _mm512_set1_epi64((long long)(c))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define AS_D(v) _mm512_castsi512_pd(v)
#define AS_I(v) _mm512_castpd_si512(v)

/* Large path (.LBL_1_12): sin(|x|) by Payne-Hanek reduction for all lanes, with the sign of x, merged into
 * res for the lanes in mask large; *nonfinite = lanes with x = +-inf or NaN. */
__attribute__((target("avx512f,avx512dq"), noinline)) static __m512d
sin_large(__m512d x, __m512d ax, __m512d res, __mmask8 large, __mmask8 *nonfinite)
{
    __m512d phi, plo, r, rlo;
    __m512i n;
    trig_reduce_avx512(ax, &phi, &plo, &n, nonfinite);
    trig_reduced_avx512(ax, phi, plo, &r, &rlo);
    const __m512d v = trig_sincos_large_tail_avx512(r, rlo, n, sin_hT, sin_hS1, sin_hS2, sin_hC1, sin_hC2, sin_hC3);
    const __m512d sign = AS_D(_mm512_xor_si512(AS_I(ax), AS_I(x)));
    return _mm512_mask_xor_pd(res, large, v, sign);
}

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_sin8_ha(__m512d x)
{
    /* y = x*16/pi + 1.5*2^52: N = round(x*16/pi) in the low bits */
    const __m512d y = _mm512_fmadd_round_pd(x, V1(sin_vInvPi), V1(sin_vRShift), RN);
    const __mmask8 mzero = _mm512_cmpeq_epi64_mask(AS_I(x), U1(sin_vMinusZero));
    const __m512d dn = _mm512_sub_round_pd(y, V1(sin_vRShift), RN);
    const __m512d ax = AS_D(_mm512_and_si512(AS_I(x), U1(sin_vAbsMask)));
    const __m512i yi = AS_I(y);
    const __m512d ch = _mm512_permutex2var_pd(_mm512_loadu_pd(sin_vCh), yi, _mm512_loadu_pd(sin_vCh + 8));
    const __m512d sh = _mm512_permutex2var_pd(_mm512_loadu_pd(sin_vSh), yi, _mm512_loadu_pd(sin_vSh + 8));
    const __m512d cl = _mm512_permutex2var_pd(_mm512_loadu_pd((const double *)sin_vClSl), yi,
                                              _mm512_loadu_pd((const double *)sin_vClSl + 8));
    /* x = -0: flip bit 4 of N so that the result is -0 */
    const __m512i ys = _mm512_mask_xor_epi64(yi, mzero, yi, U1(sin_vZeroSignBit));
    const __m512d r1 = _mm512_fnmadd_round_pd(dn, V1(sin_vPI1), x, RN);           /* x - N*PI1 */
    const __mmask8 large = _mm512_cmp_round_pd_mask(ax, V1(sin_vRangeVal), _CMP_NLE_UQ, _MM_FROUND_NO_EXC);
    const __m512i n4 = _mm512_srli_epi64(ys, 4);

    const __m512d r2 = _mm512_fnmadd_round_pd(dn, V1(sin_vPI2), r1, RN);          /* - N*PI2 */
    const __m512d r = _mm512_fnmadd_round_pd(dn, V1(sin_vPI3), r2, RN);           /* - N*PI3: R */
    __m512d e2 = _mm512_sub_round_pd(r1, r2, RN);                                 /* (N*PI2)_h */
    const __m512d rr = _mm512_mul_round_pd(r, r, RN);                             /* R^2 */
    e2 = _mm512_fnmadd_round_pd(dn, V1(sin_vPI2), e2, RN);                        /* -(N*PI2)_l */
    const __m512d e3h = _mm512_sub_round_pd(r, r2, RN);                           /* -(N*PI3)_h */
    __m512d cp = _mm512_fmadd_round_pd(V1(sin_vC4), rr, V1(sin_vC3), RN);
    const __m512d e3 = _mm512_fmadd_round_pd(V1(sin_vPI3), dn, e3h, RN);          /* (N*PI3)_l */
    const __m512d rlo = _mm512_sub_round_pd(e2, e3, RN);                          /* R_low */
    __m512d sp = _mm512_fmadd_round_pd(V1(sin_vS4), rr, V1(sin_vS3), RN);
    cp = _mm512_fmadd_round_pd(rr, cp, V1(sin_vC2), RN);
    const __m512d sl = AS_D(_mm512_slli_epi64(AS_I(cl), 32));                     /* Sl */
    cp = _mm512_fmadd_round_pd(rr, cp, V1(sin_vC1), RN);                          /* (cos R - 1)/R^2 */
    sp = _mm512_fmadd_round_pd(rr, sp, V1(sin_vS2), RN);
    const __m512d hi = _mm512_fmadd_round_pd(r, ch, sh, RN);                      /* Sh + R*Ch */
    sp = _mm512_fmadd_round_pd(rr, sp, V1(sin_vS1), RN);                          /* (sin R - R)/R^3 */
    __m512d rch_l = _mm512_sub_round_pd(hi, sh, RN);                              /* (R*Ch)_high */
    rch_l = _mm512_fmsub_round_pd(ch, r, rch_l, RN);                              /* (R*Ch)_low */
    const __m512d csr = _mm512_fnmadd_round_pd(r, sh, ch, RN);                    /* Ch - Sh*R */
    const __m512d rch = _mm512_mul_round_pd(ch, r, RN);                           /* R*Ch */
    const __m512d shc = _mm512_mul_round_pd(sh, cp, RN);                          /* Sh*CPoly */
    const __m512d lo1 = _mm512_fmadd_round_pd(r, cl, rch_l, RN);                  /* (R*Ch)_low + R*Cl */
    const __m512d lo2 = _mm512_fmadd_round_pd(rlo, csr, sl, RN);                  /* Sl + R_low*(Ch - Sh*R) */
    const __m512i sign = _mm512_slli_epi64(n4, 63);
    __m512d q = _mm512_fmadd_round_pd(rch, sp, shc, RN);                          /* Sh*CPoly + R*Ch*SPoly */
    const __m512d low = _mm512_add_round_pd(lo1, lo2, RN);
    q = _mm512_fmadd_round_pd(rr, q, low, RN);
    __m512d res = AS_D(_mm512_xor_si512(AS_I(_mm512_add_round_pd(hi, q, RN)), sign));

    if (large) {
        __mmask8 special;
        res = sin_large(x, ax, res, large, &special);
        if (special) {
            double xin[8], out[8];
            _mm512_storeu_pd(xin, x);
            _mm512_storeu_pd(out, res);
            for (int i = 0; i < 8; i++) {
                if (special & (1u << i)) {
                    svml_dsin_ha_cout_rare(&xin[i], &out[i]);
                }
            }
            res = _mm512_loadu_pd(out);
        }
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_sin_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_sin8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_sin8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
