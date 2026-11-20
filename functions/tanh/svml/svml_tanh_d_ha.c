/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision tanh, high-accuracy variant
 * (svml_z0_tanh_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), with s = sign(x), a = |x|:
 *   j  = clamp(((hi(a) & 0x7ff80000) - 0x3fc00000) >> 19, 0, 15)   (hi = high 32 bits): subinterval
 *        of a, 0 below 0.1875, then half-binades up to [16, 24), 15 from 24 on (saturation)
 *   y  = a - B[j]                   (B = midpoint of the subinterval; 0 for j = 0 and 15)
 *   tanh(a) = P0h + (P1 y + (y (y (P2 + P3 y + ... + P16 y^14)) + P0l))
 *        19 coefficient tables of 16 entries, looked up with two-source permutes; for j = 15 the
 *        polynomial is 1 + 0 y + ...
 *   result = s * tanh(a) (sign bit or-ed in)
 * Lanes with |x| >= 2^1023 * 1.5 (high word & 0x7ff80000 above 0x7fe00000; includes inf and NaN) are
 * handled by a scalar routine: +-1 for large finite x and +-inf, quiet NaN for NaN.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_tanh8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../tanh_d.h"

#if defined(__FAST_MATH__)
#error "svml_tanh_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_tanh_d_ha_constants.inc"

static inline uint64_t
dbits(double d)
{
    uint64_t u;
    memcpy(&u, &d, sizeof u);
    return u;
}

/*
 * __svml_dtanh_ha_cout_rare_internal: scalar tanh for the lanes flagged by the vector path (biased
 * exponent 2046 or 2047). Always returns 0. Statement order mirrors the assembly.
 *   NaN:                      x + x (quiet NaN)
 *   +-inf, exponent 2046:     +-1
 *   other x (never flagged):  x * (1 + x), as in the assembly
 */
int
svml_dtanh_ha_cout_rare(const double *px, double *pr)
{
    const uint64_t xb = dbits(*px);
    const double one = sOnes[xb >> 63];  /* +1 or -1 */
    const uint32_t ex = (uint32_t)(xb >> 52) & 0x7ff;

    if (ex == 0x7ff) {
        if ((xb & 0x000fffffffffffffULL) != 0) {  /* NaN */
            *pr = *px + *px;
        }
        else {
            *pr = one;
        }
        return 0;
    }
    if (ex == 0x7fe) {
        *pr = one;
        return 0;
    }
    const double x = *px;
    const double c = sOnes[0] + x;
    *pr = x * c;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define VB(b) _mm512_castsi512_pd(_mm512_set1_epi64((long long)(b)))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
/* Optimization barrier: keeps the compiler from merging the bitwise operations into different ones
 * (ternlog), so that the instruction mix stays that of the assembly (make tanh-insns). No effect on
 * results. */
#define OPAQUE(v) __asm__("" : "+v"(v))

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_tanh8_ha(__m512d x)
{
    /* subinterval index from the high words */
    const __m256i hi = _mm512_cvtepi64_epi32(_mm512_srli_epi64(_mm512_castpd_si512(x), 32));
    const __m512d ax = _mm512_and_pd(x, VB(vAbsMaskBits));
    __m512d sign = _mm512_and_pd(x, VB(vSignMaskBits));
    OPAQUE(sign);
    const __m256i hm = _mm256_and_si256(hi, _mm256_set1_epi32((int)vHiMask));
    __m256i t = _mm256_sub_epi32(hm, _mm256_set1_epi32((int)vHiBase));
    __m256 zero;
    __asm__("vxorps %0, %0, %0" : "=x"(zero));  /* the assembly's zeroing idiom (GCC would use vpxor) */
    t = _mm256_max_epi32(t, _mm256_castps_si256(zero));
    t = _mm256_min_epi32(t, _mm256_set1_epi32((int)vHiMax));
    t = _mm256_srli_epi32(t, 19);
    const __m512i j = _mm512_cvtepu32_epi64(t);
#define T(k) _mm512_permutex2var_pd(_mm512_loadu_pd(vTab[k]), j, _mm512_loadu_pd(vTab[k] + 8))

    const __m512d y = _mm512_sub_round_pd(ax, T(0), RN);  /* |x| - B */
    __m512d p = _mm512_fmadd_round_pd(y, T(18), T(17), RN);
    p = _mm512_fmadd_round_pd(y, p, T(16), RN);
    /* callout lanes: exponent field 0x7fe with the top mantissa bit set, or 0x7ff */
    const int callout = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(hm, _mm256_set1_epi32((int)vHugeHi))));
    p = _mm512_fmadd_round_pd(y, p, T(15), RN);
    p = _mm512_fmadd_round_pd(y, p, T(14), RN);
    p = _mm512_fmadd_round_pd(y, p, T(13), RN);
    p = _mm512_fmadd_round_pd(y, p, T(12), RN);
    p = _mm512_fmadd_round_pd(y, p, T(11), RN);
    p = _mm512_fmadd_round_pd(y, p, T(10), RN);
    p = _mm512_fmadd_round_pd(y, p, T(9), RN);
    p = _mm512_fmadd_round_pd(y, p, T(8), RN);
    p = _mm512_fmadd_round_pd(y, p, T(7), RN);
    p = _mm512_fmadd_round_pd(y, p, T(6), RN);
    p = _mm512_fmadd_round_pd(y, p, T(5), RN);
    p = _mm512_fmadd_round_pd(y, p, T(4), RN);       /* P2 + P3 y + ... + P16 y^14 */
    __m512d q = _mm512_mul_round_pd(p, y, RN);
    q = _mm512_fmadd_round_pd(y, q, T(2), RN);       /* y (y (...)) + P0l */
    p = _mm512_fmadd_round_pd(y, T(3), q, RN);       /* P1 y + ... */
    __m512d res = _mm512_add_round_pd(p, T(1), RN);  /* + P0h */
    res = _mm512_or_pd(res, sign);
#undef T

    if (callout) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1 << i)) {
                svml_dtanh_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_tanh_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_tanh8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 k = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, k, svml_tanh8_ha(_mm512_maskz_loadu_pd(k, x + i)));
    }
}

#endif /* __x86_64__ */
