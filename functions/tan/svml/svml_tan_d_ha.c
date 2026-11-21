/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision tan, high-accuracy variant
 * (svml_z0_tan_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Reduction, main path (|x| <= 65536.0035, 8 lanes):
 *   y  = x*16/pi + 1.5*2^52 (round to nearest): N = round(x*16/pi), j = N mod 16
 *   R + Rlo = x - N*(PI1 + PI2 + PI3)          (|R| <= pi/32; Rlo collects the rounding errors)
 * Reduction, large path (|x| > 65536.0035, +-inf, NaN; all lanes recomputed if any lane needs it): the
 *   Payne-Hanek reduction shared with sin and cos (functions/sin/svml/svml_trig_reduction.h) gives
 *   x = M*pi/256 + r + rlo (M mod 512); M = 16*k + i with k = round(M/16) (ties by the sign of r),
 *   |i| <= 8, and R + Rlo = r + i*pi/256 + rlo (pi/256 in two parts), j = k mod 16.
 * Evaluation (both paths), with Th + Tl = tan(j*pi/16) from 16-entry tables (vpermt2pd; for j = 8,
 * the pole, Th = -DBL_MAX and Tl = -2^971, which makes the quotient below -1/tan R):
 *   P + Plo = tan R = R + R^3*poly(R^2) + Rlo  (degree-4 polynomial in R^2, exact split)
 *   Num + Nlo = Th + P [+ Tl + Plo]             (two-sum)
 *   D - Dlo  = 1 - Th*P [- Th*Plo - Tl*P]       (with the exact error of the FMA)
 *   tan x = (Num + Nlo)/(D - Dlo): rcp = vrcp14pd(D) refined by one Newton-Raphson step on D - Dlo,
 *   Q = Num*rcp, then one correction step of the quotient: Q - rcp*((Q*D - Num) - Q*Dlo - Nlo).
 * Lanes with x = +-inf or NaN are then recomputed by the scalar callout (tan(+-inf) = NaN, invalid).
 *
 * vrcp14pd's result is specified by Intel's reference implementation (recip14.c, see common/include/hwymath/rcp14_ref.h);
 * the final result depends on its exact bits.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_tan8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../tan_d.h"

#if defined(__FAST_MATH__)
#error "svml_tan_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_tan_d_ha_constants.inc"

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
 * __svml_dtan_ha_cout_rare_internal: called for the lanes with x = +-inf or NaN. tan(+-inf) = 0*x =
 * NaN with the invalid flag (return code 1), tan(NaN) = x*x (quiet NaN, invalid for signaling NaNs).
 * Other inputs: returns 0 without writing *pr. The x86 result of 0*inf is the default NaN
 * 0xfff8000000000000; it is stored explicitly so that every architecture gives the same bits.
 */
int
svml_dtan_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    if (((xb >> 52) & 0x7ff) != 0x7ff) {
        return 0;
    }
    if ((uint32_t)xb == 0 && (uint32_t)(xb >> 32 & 0x7fffffff) == 0x7ff00000) {  /* +-inf */
        volatile double zero = 0.0;  /* __dtan_ha_Tab[0]; multiplied at run time: raises invalid */
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

/* Large path (.LBL_1_12): Payne-Hanek reduction of x for all lanes; replaces the table index *idx and
 * R + Rlo (*r, *rlo) in the lanes of mask large; *nonfinite = lanes with x = +-inf or NaN. */
__attribute__((target("avx512f,avx512dq"), noinline)) static void
tan_large(__m512d x, __mmask8 large, __m512i *idx, __m512d *r, __m512d *rlo, __mmask8 *nonfinite)
{
    __m512d phi, plo, z, zlo;
    __m512i n;
    trig_reduce_avx512(x, &phi, &plo, &n, nonfinite);
    trig_reduced_avx512(x, phi, plo, &z, &zlo);
    /* M = 16*k + i: k = (M + 8 - sign(z)) >> 4, i = M - 16*k (32-bit integer lanes, as in the assembly) */
    const __m256i m = _mm512_cvtepi64_epi32(n);
    const __m256i zsign = _mm256_srli_epi32(_mm512_cvtepi64_epi32(_mm512_srli_epi64(AS_I(z), 32)), 31);
    const __m256i k = _mm256_srli_epi32(_mm256_add_epi32(m, _mm256_sub_epi32(_mm256_set1_epi32((int)tan_lHalf), zsign)), 4);
    const __m512d i = _mm512_cvtepi32_pd(_mm256_sub_epi32(m, _mm256_slli_epi32(k, 4)));
    /* R + Rlo = z + i*pi/256 + zlo */
    const __m512d t = _mm512_fmadd_round_pd(V1(tan_lPi256Hi), i, z, RN);
    const __m512d u = _mm512_fnmadd_round_pd(i, V1(tan_lPi256Hi), t, RN);
    *r = _mm512_mask_blend_pd(large, *r, t);
    const __m512d w = _mm512_fmadd_round_pd(V1(tan_lPi256Lo), i, _mm512_sub_round_pd(z, u, RN), RN);
    *idx = _mm512_mask_cvtepu32_epi64(*idx, large, k);
    *rlo = _mm512_mask_add_round_pd(*rlo, large, w, zlo, RN);
}

/* Evaluation (.LBL_1_2): tan(j*pi/16 + R + Rlo), j = idx mod 16. A separate function only so that its
 * instruction mix can be compared with the corresponding part of the assembly (make tan-insns). */
__attribute__((target("avx512f,avx512dq"), noinline)) __m512d
svml_tan8_ha_eval(__m512i idx, __m512d r, __m512d rlo)
{
    const __m512d rr = _mm512_mul_round_pd(r, r, RN);
    const __m512d tl = _mm512_permutex2var_pd(_mm512_loadu_pd(tan_vTl), idx, _mm512_loadu_pd(tan_vTl + 8));
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(tan_vTh), idx, _mm512_loadu_pd(tan_vTh + 8));
    __m512d p = _mm512_fmadd_round_pd(V1(tan_vP5), rr, V1(tan_vP4), RN);
    p = _mm512_fmadd_round_pd(rr, p, V1(tan_vP3), RN);
    p = _mm512_fmadd_round_pd(rr, p, V1(tan_vP2), RN);
    p = _mm512_fmadd_round_pd(rr, p, V1(tan_vP1), RN);
    __m512d q = _mm512_mul_round_pd(p, r, RN);
    q = _mm512_fnmsub_round_pd(rr, q, rlo, RN);                    /* -(R^3*poly) - Rlo */
    const __m512d pt = _mm512_sub_round_pd(r, q, RN);              /* P = tan R, rounded */
    const __m512d plo = _mm512_sub_round_pd(_mm512_sub_round_pd(r, pt, RN), q, RN);   /* its low part */

    /* numerator Num + Nlo = Th + Tl + P + Plo */
    const __m512d num = _mm512_add_round_pd(pt, th, RN);
    __m512d nlo = _mm512_sub_round_pd(pt, _mm512_sub_round_pd(num, th, RN), RN);
    nlo = _mm512_add_round_pd(nlo, tl, RN);
    /* denominator D - Dlo = 1 - (Th + Tl)*(P + Plo) */
    const __m512d d = _mm512_fnmadd_round_pd(pt, th, V1(tan_vOne), RN);
    nlo = _mm512_add_round_pd(nlo, plo, RN);
    __m512d dlo = _mm512_sub_round_pd(d, V1(tan_vOne), RN);
    dlo = _mm512_fmadd_round_pd(pt, th, dlo, RN);
    dlo = _mm512_fmadd_round_pd(plo, th, dlo, RN);
    dlo = _mm512_fmadd_round_pd(pt, tl, dlo, RN);

    /* (Num + Nlo)/(D - Dlo) by Newton-Raphson iterations instead of a division */
    __m512d rcp = _mm512_rcp14_pd(d);
    __m512d e = _mm512_fnmadd_round_pd(d, rcp, V1(tan_vOne), RN);
    e = _mm512_fmadd_round_pd(dlo, rcp, e, RN);
    rcp = _mm512_fmadd_round_pd(e, rcp, rcp, RN);
    const __m512d quo = _mm512_mul_round_pd(rcp, num, RN);
    __m512d g = _mm512_fmsub_round_pd(quo, d, num, RN);
    g = _mm512_fnmadd_round_pd(quo, dlo, g, RN);
    g = _mm512_sub_round_pd(g, nlo, RN);
    return _mm512_fnmadd_round_pd(rcp, g, quo, RN);
}

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_tan8_ha(__m512d x)
{
    /* y = x*16/pi + 1.5*2^52: N = round(x*16/pi) in the low bits */
    const __m512d y = _mm512_fmadd_round_pd(x, V1(tan_vInvPi), V1(tan_vRShift), RN);
    const __m512d ax = AS_D(_mm512_and_si512(AS_I(x), U1(tan_vAbsMask)));
    const __m512d dn = _mm512_sub_round_pd(y, V1(tan_vRShift), RN);
    const __mmask8 large = _mm512_cmp_round_pd_mask(ax, V1(tan_vRangeVal), _CMP_NLE_UQ, _MM_FROUND_NO_EXC);
    const __m512d r1 = _mm512_fnmadd_round_pd(dn, V1(tan_vPI1), x, RN);           /* x - N*PI1 */
    const __m512d r2 = _mm512_fnmadd_round_pd(dn, V1(tan_vPI2), r1, RN);          /* - N*PI2 */
    __m512d r = _mm512_fnmadd_round_pd(dn, V1(tan_vPI3), r2, RN);                 /* - N*PI3: R */
    __m512d e2 = _mm512_sub_round_pd(r1, r2, RN);
    const __m512d e3h = _mm512_sub_round_pd(r, r2, RN);
    e2 = _mm512_fnmadd_round_pd(dn, V1(tan_vPI2), e2, RN);
    const __m512d e3 = _mm512_fmadd_round_pd(V1(tan_vPI3), dn, e3h, RN);
    __m512d rlo = _mm512_sub_round_pd(e2, e3, RN);                                /* R_low */
    __m512i idx = AS_I(y);

    __mmask8 special = 0;
    if (large) {
        tan_large(x, large, &idx, &r, &rlo, &special);
    }
    __m512d res = svml_tan8_ha_eval(idx, r, rlo);
    if (special) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (special & (1u << i)) {
                svml_dtan_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_tan_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_tan8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_tan8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
