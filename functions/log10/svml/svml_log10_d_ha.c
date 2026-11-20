/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision log10, high-accuracy variant
 * (svml_z0_log10_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), x = 2^E * m with m in [1, 2) (vgetexppd/vgetmantpd, also for denormals):
 *   rcp = vrcp14pd(m) rounded to a multiple of 1/32 (vrndscalepd): rcp = k/32, k in 16..32
 *   R   = rcp*m - 1                    (one FMA, exact)
 *   if rcp < 0.75: E += 1              (the table then holds -log10(2 rcp) instead of -log10(rcp))
 *   j   = top 4 mantissa bits of rcp
 *   log10(x) = E*log10(2) + T[j] + log10(1+R),  T[j] = Th[j] + Tl[j]
 *   log10(1+R) = c1h*R + (c1l*R + c2 R^2 + ... + c9 R^9),  c1h + c1l = 1/ln10
 *   KT = E*L2H + Th[j] (exact),  H = KT + c1h*R (FMA),  (c1h*R)_h = H - KT,  (c1h*R)_l = c1h*R - (c1h*R)_h
 *   result = H + (R*p + ((E*L2L + Tl[j]) + (c1h*R)_l))
 * Lanes with x = +-0, +-inf or x < 0 (vfpclasspd 0x5e; NaN and denormals stay in the vector path) are
 * recomputed by a scalar routine, which also returns the error code and raises the FP flags.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_log108_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../log10_d.h"

#if defined(__FAST_MATH__)
#error "svml_log10_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_log10_d_ha_constants.inc"

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
 * __svml_dlog10_ha_cout_rare_internal: scalar log with a 65-entry table {rcp, -log(rcp) hi, lo} indexed by
 * the top 6 mantissa bits (rounded), a degree-8 polynomial, and a separate polynomial for |x-1| <= 2^-6.6.
 * Statement order mirrors the assembly.
 */
int
svml_dlog10_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);

    if (((xb >> 48) & 0x7ff0) == 0x7ff0) {  /* inf or NaN */
        if ((xb >> 63) && (xb & 0x000fffffffffffffULL) == 0) {  /* -inf: NaN, invalid */
            volatile double z = sZero;
            *pr = copysign(z / z, -1.0);  /* x86's default NaN (negative) on every architecture */
            return 1;
        }
        *pr = x * x;  /* +inf -> +inf, NaN -> quiet NaN */
        return 0;
    }
    double x2 = x;
    int ecorr = 0;
    if (((dbits(x2) >> 48) & 0x7ff0) == 0) {  /* zero or denormal: scale by 2^60 */
        x2 = x2 * sTwo60;
        ecorr = -60;
    }
    if (!(x2 > sZero)) {  /* x <= 0 */
        if (x2 == sZero) {  /* -inf, divide-by-zero */
            volatile double z = sZero;
            *pr = -sOne / z;
            return 2;
        }
        volatile double z = sZero;  /* x < 0: NaN, invalid */
        *pr = copysign(z / z, -1.0);  /* x86's default NaN (negative) on every architecture */
        return 1;
    }
    const double t = x2 - sPacketOne;
    if (!(fabs(t) > sNear1)) {  /* x close to 1: s = t*sScale, log10(1+t) = s + s*(p0 + ... + p7 s^7) */
        const double ts = t * sScale;
        double p = sP7 * ts;
        p = p + sP6;
        p = p * ts;
        p = p + sP5;
        p = p * ts;
        p = p + sP4;
        p = p * ts;
        p = p + sP3;
        p = p * ts;
        p = p + sP2;
        p = p * ts;
        p = p + sP1;
        p = p * ts;
        p = p + sP0;
        p = p * ts;
        *pr = p + ts;
        return 0;
    }

    /* y = mantissa of x2 in [1, 2), e = its exponent */
    const uint64_t x2b = dbits(x2);
    const double y = bitsd((x2b & 0x800fffffffffffffULL) | 0x3ff0000000000000ULL);
    const int e = (int)((x2b >> 52) & 0x7ff) - 1023 + ecorr;
    const double a = y + sIdxShifter;  /* low bits: index of the nearest rcp */
    const uint32_t idx = (uint32_t)dbits(a) & 127;
    double yh = y + sRndShifter;  /* y rounded to 29 fractional bits */
    yh = yh - sRndShifter;
    const double *T = &sTab[3 * idx];  /* T[0] = rcp*sScale */
    const double ed = (double)e;
    const double yl = y - yh;
    double u = yh * T[0];
    double hi = sL2H * ed;
    u = u - sScale;
    const double v = yl * T[0];
    double lo = ed * sL2L;
    hi = hi + T[1];
    lo = lo + T[2];
    hi = hi + u;
    const double r = u + v;
    double p = sP7 * r;
    p = p + sP6;
    p = p * r;
    p = p + sP5;
    p = p * r;
    p = p + sP4;
    p = p * r;
    p = p + sP3;
    p = p * r;
    p = p + sP2;
    p = p * r;
    p = p + sP1;
    p = p * r;
    p = p + sP0;
    const double pu = u * p;
    p = p * v;
    lo = lo + p;
    lo = lo + v;
    *pr = hi + (pu + lo);
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_log108_ha(__m512d x)
{
    const __m512d m = _mm512_getmant_round_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_nan, SAE);  /* [1,2), NaN for x < 0 */
    __m512d e = _mm512_getexp_round_pd(x, SAE);
    const __m512d rcp0 = _mm512_rcp14_pd(m);
    const __mmask8 special = _mm512_fpclass_pd_mask(x, 0x5e);  /* +-0, +-inf, negative */
    const __m512d rcp = _mm512_roundscale_round_pd(rcp0, 0x58, SAE);  /* RN to 5 fractional bits */
    const __m512d r = _mm512_fmsub_round_pd(rcp, m, V1(vOne), RN);   /* R = rcp*m - 1 */
    const __mmask8 lt34 = _mm512_cmp_round_pd_mask(rcp, V1(vThree4), _CMP_LT_OQ, SAE);
    const __m512d p67 = _mm512_fmadd_round_pd(V1(vC7), r, V1(vC6), RN);
    const __m512d p89 = _mm512_fmadd_round_pd(V1(vC9), r, V1(vC8), RN);
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(rcp), 48);
    const __m512d r4 = _mm512_mul_round_pd(r2, r2, RN);
    const __m512d p69 = _mm512_fmadd_round_pd(r2, p89, p67, RN);
    e = _mm512_mask_add_round_pd(e, lt34, e, V1(vOne), RN);          /* E + 1 if rcp < 0.75 */
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));
    const __m512d tl = _mm512_permutex2var_pd(_mm512_loadu_pd(vTl), idx, _mm512_loadu_pd(vTl + 8));
    const __m512d kt = _mm512_fmadd_round_pd(e, V1(vL2H), th, RN);   /* E*log10(2)_high + Th */
    const __m512d ktl = _mm512_fmadd_round_pd(V1(vL2L), e, tl, RN);  /* E*log10(2)_low + Tl */
    const __m512d p23 = _mm512_fmadd_round_pd(V1(vC3), r, V1(vC2), RN);
    const __m512d p45 = _mm512_fmadd_round_pd(V1(vC5), r, V1(vC4), RN);
    const __m512d h = _mm512_fmadd_round_pd(r, V1(vC1H), kt, RN);    /* K + Th + (c1h*R)_h */
    const __m512d p25 = _mm512_fmadd_round_pd(r2, p45, p23, RN);
    __m512d crh = _mm512_sub_round_pd(h, kt, RN);                     /* (c1h*R)_h */
    __m512d p = _mm512_fmadd_round_pd(r4, p69, p25, RN);              /* polynomial */
    crh = _mm512_mask_xor_pd(crh, special, crh, crh);
    p = _mm512_fmadd_round_pd(r, p, V1(vC1L), RN);
    const __m512d crl = _mm512_fmsub_round_pd(r, V1(vC1H), crh, RN);  /* (c1h*R)_l */
    const __m512d lo = _mm512_add_round_pd(ktl, crl, RN);
    p = _mm512_fmadd_round_pd(r, p, lo, RN);
    __m512d res = _mm512_add_round_pd(p, h, RN);

    if (special) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (special & (1u << i)) {
                svml_dlog10_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_log10_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_log108_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        /* padding lanes are loaded as 1.0 (not 0, which would take the scalar callout) */
        _mm512_mask_storeu_pd(y + i, m, svml_log108_ha(_mm512_mask_loadu_pd(_mm512_set1_pd(1.0), m, x + i)));
    }
}

#endif /* __x86_64__ */
