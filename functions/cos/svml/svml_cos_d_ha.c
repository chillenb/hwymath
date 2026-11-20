/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision cos, high-accuracy variant
 * (svml_z0_cos_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * The algorithm is that of sin (functions/sin/svml/svml_sin_d_ha.c), with the same tables:
 * Main path (|x| <= 2^24), 8 lanes:
 *   y  = x*16/pi + 1.5*2^52 (round to nearest): N = round(x*16/pi), j = N mod 16
 *   R + Rlo = x - N*(PI1 + PI2 + PI3)
 *   B = N*pi/16: Sh + Sl = sin B, Ch + Cl = cos B (16-entry tables, vpermt2pd)
 *   cos(B + R) = Ch - R*Sh [split exactly] + R^2*(Ch*Cp(R^2) - R*Sh*Sp(R^2))
 *                + Cl - Rlo*(Sh + Ch*R) - [(error of R*Sh) + R*Sl]
 *   bit 4 of N flips the sign.
 * Large path (|x| > 2^24, +-inf, NaN): Payne-Hanek reduction of x (with its sign) modulo pi/256 and
 * the 512-entry table of cos/sin(M*pi/256), shared with sin (functions/sin/svml/svml_trig_reduction.h).
 * Lanes with x = +-inf or NaN are then recomputed by the scalar callout (cos(+-inf) = NaN, invalid).
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_cos8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../cos_d.h"

#if defined(__FAST_MATH__)
#error "svml_cos_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_cos_d_ha_constants.inc"

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
 * __svml_dcos_ha_cout_rare_internal: called for the lanes with x = +-inf or NaN. cos(+-inf) = 0*x =
 * NaN with the invalid flag (return code 1), cos(NaN) = x*x (quiet NaN, invalid for signaling NaNs).
 * Other inputs: returns 0 without writing *pr. The x86 result of 0*inf is the default NaN
 * 0xfff8000000000000; it is stored explicitly so that every architecture gives the same bits.
 */
int
svml_dcos_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    if (((xb >> 52) & 0x7ff) != 0x7ff) {
        return 0;
    }
    if ((uint32_t)xb == 0 && (uint32_t)(xb >> 32 & 0x7fffffff) == 0x7ff00000) {  /* +-inf */
        volatile double zero = 0.0;  /* __dcos_ha_CoutTab[0]; multiplied at run time: raises invalid */
        volatile double r = x * zero;
        (void)r;
        *pr = bitsd(0xfff8000000000000ULL);
        return 1;
    }
    *pr = x * x;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#include "functions/sin/svml/svml_trig_reduction.h"

#define V1(c) _mm512_set1_pd(c)
#define U1(c) _mm512_set1_epi64((long long)(c))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define AS_D(v) _mm512_castsi512_pd(v)
#define AS_I(v) _mm512_castpd_si512(v)

/* Large path (.LBL_1_12): cos(x) by Payne-Hanek reduction for all lanes, merged into res for the lanes
 * in mask large; *nonfinite = lanes with x = +-inf or NaN. */
__attribute__((target("avx512f,avx512dq"), noinline)) static __m512d
cos_large(__m512d x, __m512d res, __mmask8 large, __mmask8 *nonfinite)
{
    __m512d phi, plo, r, rlo;
    __m512i n;
    trig_reduce_avx512(x, &phi, &plo, &n, nonfinite);
    trig_reduced_avx512(x, phi, plo, &r, &rlo);
    const __m512d v = trig_sincos_large_tail_avx512(r, rlo, n, cos_hT, cos_hS1, cos_hS2, cos_hC1, cos_hC2, cos_hC3);
    return _mm512_mask_mov_pd(res, large, v);
}

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_cos8_ha(__m512d x)
{
    /* y = x*16/pi + 1.5*2^52: N = round(x*16/pi) in the low bits */
    const __m512d y = _mm512_fmadd_round_pd(x, V1(cos_vInvPi), V1(cos_vRShift), RN);
    const __m512d ax = AS_D(_mm512_and_si512(AS_I(x), U1(cos_vAbsMask)));
    const __m512d dn = _mm512_sub_round_pd(y, V1(cos_vRShift), RN);
    const __m512i yi = AS_I(y);
    const __m512d sh = _mm512_permutex2var_pd(_mm512_loadu_pd(cos_vSh), yi, _mm512_loadu_pd(cos_vSh + 8));
    const __m512d ch = _mm512_permutex2var_pd(_mm512_loadu_pd(cos_vCh), yi, _mm512_loadu_pd(cos_vCh + 8));
    const __m512d cl = _mm512_permutex2var_pd(_mm512_loadu_pd((const double *)cos_vClSl), yi,
                                              _mm512_loadu_pd((const double *)cos_vClSl + 8));
    const __m512d r1 = _mm512_fnmadd_round_pd(dn, V1(cos_vPI1), x, RN);           /* x - N*PI1 */
    const __mmask8 large = _mm512_cmp_round_pd_mask(ax, V1(cos_vRangeVal), _CMP_NLE_UQ, _MM_FROUND_NO_EXC);
    const __m512i n4 = _mm512_srli_epi64(yi, 4);

    const __m512d r2 = _mm512_fnmadd_round_pd(dn, V1(cos_vPI2), r1, RN);          /* - N*PI2 */
    const __m512d r = _mm512_fnmadd_round_pd(dn, V1(cos_vPI3), r2, RN);           /* - N*PI3: R */
    __m512d e2 = _mm512_sub_round_pd(r1, r2, RN);                                 /* (N*PI2)_h */
    const __m512d rr = _mm512_mul_round_pd(r, r, RN);                             /* R^2 */
    e2 = _mm512_fnmadd_round_pd(dn, V1(cos_vPI2), e2, RN);                        /* -(N*PI2)_l */
    const __m512d e3h = _mm512_sub_round_pd(r, r2, RN);                           /* -(N*PI3)_h */
    const __m512d e3 = _mm512_fmadd_round_pd(V1(cos_vPI3), dn, e3h, RN);          /* (N*PI3)_l */
    __m512d cp = _mm512_fmadd_round_pd(V1(cos_vC4), rr, V1(cos_vC3), RN);
    const __m512d rlo = _mm512_sub_round_pd(e2, e3, RN);                          /* R_low */
    __m512d sp = _mm512_fmadd_round_pd(V1(cos_vS4), rr, V1(cos_vS3), RN);
    cp = _mm512_fmadd_round_pd(rr, cp, V1(cos_vC2), RN);
    const __m512d sl = AS_D(_mm512_slli_epi64(AS_I(cl), 32));                     /* Sl */
    cp = _mm512_fmadd_round_pd(rr, cp, V1(cos_vC1), RN);                          /* (cos R - 1)/R^2 */
    sp = _mm512_fmadd_round_pd(rr, sp, V1(cos_vS2), RN);
    const __m512d hi = _mm512_fnmadd_round_pd(r, sh, ch, RN);                     /* Ch - R*Sh */
    sp = _mm512_fmadd_round_pd(rr, sp, V1(cos_vS1), RN);                          /* (sin R - R)/R^3 */
    __m512d rsh_l = _mm512_sub_round_pd(ch, hi, RN);                              /* (R*Sh)_high */
    rsh_l = _mm512_fmsub_round_pd(r, sh, rsh_l, RN);                              /* (R*Sh)_low */
    const __m512d sc = _mm512_fmadd_round_pd(r, ch, sh, RN);                      /* Sh + Ch*R */
    const __m512d rsh = _mm512_mul_round_pd(sh, r, RN);                           /* R*Sh */
    const __m512d chc = _mm512_mul_round_pd(ch, cp, RN);                          /* Ch*CPoly */
    const __m512d lo1 = _mm512_fmadd_round_pd(r, sl, rsh_l, RN);                  /* (R*Sh)_low + R*Sl */
    const __m512d lo2 = _mm512_fnmadd_round_pd(rlo, sc, cl, RN);                  /* Cl - R_low*(Sh + Ch*R) */
    const __m512i sign = _mm512_slli_epi64(n4, 63);
    __m512d q = _mm512_fnmadd_round_pd(rsh, sp, chc, RN);                         /* Ch*CPoly - R*Sh*SPoly */
    const __m512d low = _mm512_sub_round_pd(lo2, lo1, RN);
    q = _mm512_fmadd_round_pd(rr, q, low, RN);
    __m512d res = AS_D(_mm512_xor_si512(AS_I(_mm512_add_round_pd(hi, q, RN)), sign));

    if (large) {
        __mmask8 special;
        res = cos_large(x, res, large, &special);
        if (special) {
            double xin[8], out[8];
            _mm512_storeu_pd(xin, x);
            _mm512_storeu_pd(out, res);
            for (int i = 0; i < 8; i++) {
                if (special & (1u << i)) {
                    svml_dcos_ha_cout_rare(&xin[i], &out[i]);
                }
            }
            res = _mm512_loadu_pd(out);
        }
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_cos_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_cos8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_cos8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
