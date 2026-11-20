/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision atan, high-accuracy variant
 * (svml_z0_atan_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results are bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes, no scalar callout), with a = |x| and the sign of x restored at the end:
 *   a < 7.875:  x0 = j/4 = RN(a + 1.5*2^50) - 1.5*2^50 (a rounded to a multiple of 1/4, j = 0..31 from the
 *               low mantissa bits), Diff = a - x0 (vreducepd, exact), Y = 1 + x0*a (one FMA),
 *               atan(a) = atan(j/4) + atan(Diff/Y), atan(j/4) = AtanHi[j] + AtanLo[j] (two 16-entry
 *               tables each, vpermt2pd, selected by j >= 16);
 *   a >= 7.875: Diff = -1, Y = min(a, 2^128), atan(a) = pi/2 + atan(-1/Y).
 *   R + R_low = Diff/Y: Rcp = rcp14(Y) (vrcp14pd, |rel err| < 2^-14), eps = 1 - Y*Rcp,
 *               Rcp2 = Rcp(1 + eps)(1 + eps^2) (two FMAs), R = RN(Diff*Rcp2),
 *               R_low = (Diff*Rcp2 - R) + R*(1 - Y*Rcp2) [- (x0*a - (Y - 1))*Rcp2*R for a < 7.875: the
 *               rounding error of Y];
 *   atan(a) = (AtanHi + R) + (R^3 P(R^2) + (R_low + AtanLo + err(AtanHi + R))),  P: degree 5.
 * NaN inputs give the quieted NaN; +-inf give +-pi/2 (rounded).
 *
 * The scalar routine (__svml_datan_ha_cout_rare_internal, a complete table-driven scalar atan) is not
 * called by this vector code; it is translated mechanically, instruction by instruction, by
 * functions/asin/tools/asm2c.py into svml_datan_ha_cout_rare.inc and tested on its own.
 *
 * vrcp14pd results are implementation-defined (Intel's are a function of the top 16 mantissa bits, with
 * at most 16 significant bits); the result of the algorithm depends on them (Rcp2 keeps the rounding
 * error of Rcp(1 + eps)). svml_atan8_ha_portable() replaces the instruction by the deterministic
 * approximation used by the Highway version, trunc16(RN(1/Y)) (16 significant bits, relative error
 * < 2^-16 + 2^-53), as a reference for the Highway targets.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only the vector functions need AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../atan_d.h"
#include "asm_model.h"

#if defined(__FAST_MATH__)
#error "svml_atan_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_atan_d_ha_constants.inc"

/* __svml_datan_ha_cout_rare_internal (not called by the vector path) */
#include "svml_datan_ha_cout_rare.inc"

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC
#define AVX512 __attribute__((target("avx512f,avx512dq")))
/* compiler barrier without instructions: keeps the instruction mix comparable with the assembly */
#define OPAQUE(v) __asm__("" : "+v"(v))

/* trunc16(RN(1/y)): the portable replacement of vrcp14pd (see above) */
AVX512 static inline __m512d
rcp_portable(__m512d y)
{
    const __m512d r = _mm512_div_round_pd(V1(1.0), y, RN);
    const __m512d t = _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(r), _mm512_set1_epi64((long long)0xfffffff000000000ULL)));
    /* NaN: the quieted input, as from the instruction (keeps the NaN payloads of the results unchanged) */
    return _mm512_mask_add_round_pd(t, _mm512_cmp_pd_mask(y, y, _CMP_UNORD_Q), y, y, RN);
}

AVX512 static inline __attribute__((always_inline)) __m512d
atan8(__m512d x, int portable)
{
    const __m512d one = V1(vOne);
    __m512d ax = _mm512_castsi512_pd(
            _mm512_and_si512(_mm512_castpd_si512(x), _mm512_set1_epi64((long long)vAbsMaskBits)));
    OPAQUE(ax);
    const __mmask8 small = _mm512_cmp_round_pd_mask(ax, V1(vTabMax), _CMP_LT_OQ, SAE);  /* a < 7.875 */
    __m512d y = _mm512_min_round_pd(V1(vYMax), ax, SAE);
    /* Diff = a - round(4a)/4 (vreducepd imm 0x28: 2 fraction bits, round to nearest even, no precision
     * exception), or -1 */
    const __m512d diff = _mm512_mask_reduce_round_pd(V1(vMOne), small, ax, 0x28, SAE);
    const __m512d t = _mm512_add_round_pd(ax, V1(vShifter), RN);  /* low bits: j */
    __m512d sgn = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(ax), _mm512_castpd_si512(x)));
    OPAQUE(sgn);
    const __mmask8 large = (__mmask8)~small;
    const __m512d x0 = _mm512_sub_round_pd(t, V1(vShifter), RN);  /* j/4 */
    const __m512i idx = _mm512_castpd_si512(t);
    const __m512d th0 = _mm512_permutex2var_pd(_mm512_loadu_pd(vAtanHi), idx, _mm512_loadu_pd(vAtanHi + 8));
    y = _mm512_mask_fmadd_round_pd(y, small, x0, one, RN);  /* Y = a*x0 + 1 for a < 7.875 */
    const __mmask8 hi16 = _mm512_cmp_round_pd_mask(t, V1(vShifter16), _CMP_GE_OQ, SAE);  /* j >= 16 */
    __m512d rcp = portable ? rcp_portable(y) : _mm512_rcp14_pd(y);
    const __m512d th1 = _mm512_permutex2var_pd(_mm512_loadu_pd(vAtanHi + 16), idx, _mm512_loadu_pd(vAtanHi + 24));
    const __m512d ym1 = _mm512_sub_round_pd(y, one, RN);             /* (x*x0)_high */
    const __m512d xx0lo = _mm512_fmsub_round_pd(x0, ax, ym1, RN);   /* (x*x0)_low */
    const __m512d tl0 = _mm512_permutex2var_pd(_mm512_loadu_pd(vAtanLo), idx, _mm512_loadu_pd(vAtanLo + 8));
    const __m512d tl1 = _mm512_permutex2var_pd(_mm512_loadu_pd(vAtanLo + 16), idx, _mm512_loadu_pd(vAtanLo + 24));
    __m512d thi = _mm512_mask_blend_pd(hi16, th0, th1);

    /* Rcp2 = Rcp (1 + eps)(1 + eps^2), eps = 1 - Y Rcp */
    const __m512d eps = _mm512_fnmadd_round_pd(rcp, y, one, RN);
    const __m512d eps2 = _mm512_mul_round_pd(eps, eps, RN);
    rcp = _mm512_fmadd_round_pd(eps, rcp, rcp, RN);
    __m512d tlo = _mm512_mask_blend_pd(hi16, tl0, tl1);
    rcp = _mm512_fmadd_round_pd(eps2, rcp, rcp, RN);

    /* R + R_low = Diff/Y */
    const __m512d r = _mm512_mul_round_pd(rcp, diff, RN);
    const __m512d xx0lor = _mm512_mul_round_pd(xx0lo, rcp, RN);
    const __m512d eps1 = _mm512_fnmadd_round_pd(rcp, y, one, RN);   /* 1 - Y Rcp2 */
    const __m512d rlo0 = _mm512_fmsub_round_pd(rcp, diff, r, RN);   /* (Diff*Rcp2)_low */
    __m512d rlo = _mm512_fmadd_round_pd(r, eps1, rlo0, RN);
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    rlo = _mm512_mask3_fnmadd_round_pd(xx0lor, r, rlo, small, RN);  /* - (x*x0)_low*Rcp2*R */

    /* pi/2 for large a */
    thi = _mm512_mask_blend_pd(large, thi, V1(vPio2Hi));
    tlo = _mm512_mask_blend_pd(large, tlo, V1(vPio2Lo));
    const __m512d r4 = _mm512_mul_round_pd(r2, r2, RN);
    const __m512d r3 = _mm512_mul_round_pd(r2, r, RN);
    __m512d p45 = _mm512_fmadd_round_pd(V1(vA5), r2, V1(vA4), RN);
    const __m512d p23 = _mm512_fmadd_round_pd(V1(vA3), r2, V1(vA2), RN);
    const __m512d low = _mm512_add_round_pd(rlo, tlo, RN);
    const __m512d hr = _mm512_add_round_pd(thi, r, RN);
    p45 = _mm512_fmadd_round_pd(r4, p45, p23, RN);
    const __m512d hr_t = _mm512_sub_round_pd(hr, thi, RN);
    const __m512d p01 = _mm512_fmadd_round_pd(V1(vA1), r2, V1(vA0), RN);
    const __m512d hr_err = _mm512_sub_round_pd(r, hr_t, RN);        /* R - ((AtanHi + R) - AtanHi) */
    __m512d p = _mm512_fmadd_round_pd(r4, p45, p01, RN);
    const __m512d low2 = _mm512_add_round_pd(low, hr_err, RN);
    p = _mm512_fmadd_round_pd(r3, p, low2, RN);
    __m512d res = _mm512_add_round_pd(p, hr, RN);
    res = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(res), _mm512_castpd_si512(sgn)));
    return res;
}

AVX512 __m512d
svml_atan8_ha(__m512d x)
{
    return atan8(x, 0);
}

AVX512 __m512d
svml_atan8_ha_portable(__m512d x)
{
    return atan8(x, 1);
}

#define ARRAY(NAME, KERNEL)                                                                         \
    AVX512 void NAME(const double *x, double *y, size_t len)                                        \
    {                                                                                               \
        size_t i = 0;                                                                               \
        for (; i + 8 <= len; i += 8) {                                                              \
            _mm512_storeu_pd(y + i, KERNEL(_mm512_loadu_pd(x + i)));                                \
        }                                                                                           \
        if (i < len) {                                                                              \
            const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);                                   \
            _mm512_mask_storeu_pd(y + i, m, KERNEL(_mm512_maskz_loadu_pd(m, x + i)));               \
        }                                                                                           \
    }
ARRAY(svml_atan_ha_array, svml_atan8_ha)
ARRAY(svml_atan_ha_portable_array, svml_atan8_ha_portable)

#endif /* __x86_64__ */
