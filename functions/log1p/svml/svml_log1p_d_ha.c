/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision log1p, high-accuracy variant
 * (svml_z0_log1p_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes): s = 1 + x is computed as a double-double Xh + Xl (Xh = RN(1+x),
 * Xl = min(1,x) - (Xh - max(1,x)), exact), then log(Xh + Xl) like SVML's log:
 *   Xh = 2^E * m with m in [1, 2) (vgetexppd/vgetmantpd)
 *   rcp = vrcp14pd(m) rounded to a multiple of 1/32 (vrndscalepd): rcp = k/32, k in 16..32
 *   R   = rcp*m - 1 (exact) + rcp*Xl*2^-E (vscalefpd),  as R' = R + rcpXl with its error Rl
 *   if rcp < 0.75: E += 1;   j = top 4 mantissa bits of rcp
 *   log1p(x) = E*ln2 + T[j] + log(1+R'),  T[j] = Th[j] + Tl[j] = -log(rcp) (or -log(2 rcp))
 *   High = E*L2H + Th[j] (exact),  HR = High + R',  Rh = HR - High,  Rl2 = R' - Rh
 *   result = (HR + ((R'^2*p + Rl2) + (E*L2L + (Tl[j] + Rl)))) | sign(x)   (log1p(-0) = -0)
 * Lanes where Rl2 is NaN (x < -1, x = -1, +-inf, NaN) are recomputed by a scalar routine (SVML's
 * scalar log applied to 1 + x), which also returns the error code and raises the FP flags.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_log1p8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../log1p_d.h"

#if defined(__FAST_MATH__)
#error "svml_log1p_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_log1p_d_ha_constants.inc"

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
 * __svml_dlog1p_ha_cout_rare_internal: SVML's scalar log applied to x2 = RN(1 + x), with a 65-entry table {rcp, -log(rcp) hi, lo} indexed by
 * the top 6 mantissa bits (rounded), a degree-8 polynomial, and a separate polynomial for |x-1| <= 2^-6.6.
 * Statement order mirrors the assembly.
 */
int
svml_dlog1p_ha_cout_rare(const double *px, double *pr)
{
    double x2 = *px + sPacketOne;  /* 1 + x */
    const uint64_t xb = dbits(x2);

    if (((xb >> 48) & 0x7ff0) == 0x7ff0) {  /* 1 + x is inf or NaN */
        if ((xb >> 63) && (xb & 0x000fffffffffffffULL) == 0) {  /* -inf: NaN, invalid */
            volatile double z = sZero;
            *pr = copysign(z / z, -1.0);  /* x86's default NaN (negative) on every architecture */
            return 1;
        }
        *pr = x2 * x2;  /* +inf -> +inf, NaN -> quiet NaN */
        return 0;
    }
    int ecorr = 0;
    if (((dbits(x2) >> 48) & 0x7ff0) == 0) {  /* zero or denormal: scale by 2^60 */
        x2 = x2 * sTwo60;
        ecorr = -60;
    }
    if (!(x2 > sZero)) {  /* 1 + x <= 0 */
        if (x2 == sZero) {  /* -inf, divide-by-zero */
            volatile double z = sZero;
            *pr = -sOne / z;
            return 2;
        }
        volatile double z = sZero;  /* 1 + x < 0: NaN, invalid */
        *pr = copysign(z / z, -1.0);  /* x86's default NaN (negative) on every architecture */
        return 1;
    }
    const double t = x2 - sPacketOne;
    if (!(fabs(t) > sNear1)) {  /* 1 + x close to 1: log(1+t) = t + t^2*(c2 + ... + c8 t^6) */
        double p = sP6 * t;
        p = p + sP5;
        p = p * t;
        p = p + sP4;
        p = p * t;
        p = p + sP3;
        p = p * t;
        p = p + sP2;
        p = p * t;
        p = p + sP1;
        p = p * t;
        p = p + sP0;
        p = p * t;
        p = p * t;
        *pr = p + t;
        return 0;
    }

    /* y = mantissa of x2 in [1, 2), e = its exponent */
    const uint64_t x2b = dbits(x2);
    const double y = bitsd((x2b & 0x800fffffffffffffULL) | 0x3ff0000000000000ULL);
    const int e = (int)((x2b >> 52) & 0x7ff) - 1023 + ecorr;
    const double a = y + sIdxShifter;  /* low bits: index of the nearest rcp */
    const uint32_t idx = (uint32_t)dbits(a) & 127;
    double yh = y + sRndShifter;  /* y rounded to 31 fractional bits */
    yh = yh - sRndShifter;
    const double *T = &sTab[3 * idx];
    const double ed = (double)e;
    const double yl = y - yh;
    double u = yh * T[0];
    double hi = sL2H * ed;
    u = u - sPacketOne;  /* exact */
    const double v = yl * T[0];
    double lo = ed * sL2L;
    hi = hi + T[1];
    lo = lo + T[2];
    hi = hi + u;
    lo = lo + v;
    const double r = u + v;
    double p = sP6 * r;
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
    p = p * r;
    p = p * r;
    p = p + lo;
    *pr = hi + p;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC

__attribute__((target("avx512f"))) __m512d
svml_log1p8_ha(__m512d x)
{
    /* 1 + x as Xh + Xl */
    const __m512d mx = _mm512_max_round_pd(V1(vOne), x, SAE);
    const __m512d mn = _mm512_min_round_pd(V1(vOne), x, SAE);
    const __m512d xh = _mm512_add_round_pd(V1(vOne), x, RN);
    const __m512i sgn = _mm512_and_si512(_mm512_castpd_si512(x), _mm512_set1_epi64((long long)vSignMask));
    const __m512d m = _mm512_getmant_round_pd(xh, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_nan, SAE);
    __m512d e = _mm512_getexp_round_pd(xh, SAE);
    const __m512d d = _mm512_sub_round_pd(xh, mx, RN);
    const __m512d rcp0 = _mm512_rcp14_pd(m);
    const __m512d xl = _mm512_sub_round_pd(mn, d, RN);                /* Xl */
    const __m512d nege = _mm512_castsi512_pd(
            _mm512_xor_si512(_mm512_castpd_si512(e), _mm512_set1_epi64((long long)vSignMask)));
    const __m512d rcp = _mm512_roundscale_round_pd(rcp0, 0x58, SAE);  /* RN to 5 fractional bits */
    const __m512d xls = _mm512_scalef_round_pd(xl, nege, RN);         /* Xl*2^(-E) */
    const __mmask8 lt34 = _mm512_cmp_round_pd_mask(rcp, V1(vThree4), _CMP_LT_OQ, SAE);
    const __m512d r = _mm512_fmsub_round_pd(rcp, m, V1(vOne), RN);   /* R = rcp*m - 1 */
    const __m512d rxl = _mm512_mul_round_pd(rcp, xls, RN);
    e = _mm512_mask_add_round_pd(e, lt34, e, V1(vOne), RN);          /* E + 1 if rcp < 0.75 */
    const __m512d rr = _mm512_add_round_pd(r, rxl, RN);              /* R' */
    const __m512d r2 = _mm512_mul_round_pd(rr, rr, RN);
    const __m512d p89 = _mm512_fmadd_round_pd(V1(vC9), rr, V1(vC8), RN);
    __m512d rl = _mm512_sub_round_pd(rr, r, RN);
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(rcp), 48);
    rl = _mm512_sub_round_pd(rxl, rl, RN);                            /* Rl = rcpXl - (R' - R) */
    const __m512d tl = _mm512_permutex2var_pd(_mm512_loadu_pd(vTl), idx, _mm512_loadu_pd(vTl + 8));
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));
    const __m512d high = _mm512_fmadd_round_pd(V1(vL2H), e, th, RN); /* E*log(2)_high + Th */
    const __m512d tlrl = _mm512_add_round_pd(tl, rl, RN);             /* Rl + Tl */
    const __m512d p67 = _mm512_fmadd_round_pd(V1(vC7), rr, V1(vC6), RN);
    const __m512d p23 = _mm512_fmadd_round_pd(V1(vC3), rr, V1(vC2), RN);
    const __m512d hr = _mm512_add_round_pd(high, rr, RN);             /* E*log(2)_high + Th + R' */
    const __m512d p69 = _mm512_fmadd_round_pd(r2, p89, p67, RN);
    const __m512d p45 = _mm512_fmadd_round_pd(V1(vC5), rr, V1(vC4), RN);
    const __m512d r4 = _mm512_mul_round_pd(r2, r2, RN);
    const __m512d p25 = _mm512_fmadd_round_pd(r2, p45, p23, RN);
    const __m512d rh = _mm512_sub_round_pd(hr, high, RN);
    const __m512d p = _mm512_fmadd_round_pd(r4, p69, p25, RN);        /* polynomial */
    const __m512d rl2 = _mm512_sub_round_pd(rr, rh, RN);
    const __m512d q = _mm512_fmadd_round_pd(r2, p, rl2, RN);
    const __mmask8 special = _mm512_cmp_round_pd_mask(rl2, rl2, _CMP_NEQ_UQ, SAE);  /* NaN */
    const __m512d lo = _mm512_fmadd_round_pd(V1(vL2L), e, tlrl, RN); /* E*log(2)_low + Tl + Rl */
    __m512d res = _mm512_add_round_pd(q, lo, RN);
    res = _mm512_add_round_pd(hr, res, RN);
    res = _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(res), sgn));  /* log1p(-0) = -0 */

    if (special) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (special & (1u << i)) {
                svml_dlog1p_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f"))) void
svml_log1p_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_log1p8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_log1p8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
