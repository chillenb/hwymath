/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision log2, high-accuracy variant
 * (svml_z0_log2_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), x = 2^E * m with m in [1, 2) (vgetexppd/vgetmantpd, also for denormals):
 *   rcp = vrcp14pd(m) rounded to a multiple of 1/32 (vrndscalepd): rcp = k/32, k in 16..32
 *   R   = rcp*m - 1                    (one FMA, exact)
 *   if rcp < 0.75: E += 1              (the table then holds -log2(2 rcp) instead of -log2(rcp))
 *   j   = top 4 mantissa bits of rcp
 *   log2(x) = E + T[j] + log2(1+R),    T[j] = Th[j] + Tl[j]
 *   log2(1+R) = c1h*R + R*(c1l + c2 R + ... + c10 R^9),  c1h + c1l = 1/ln2
 *   KT = E + Th[j],  H = KT + c1h*R (FMA),  (c1h*R)_h = H - KT,  (c1h*R)_l = c1h*R - (c1h*R)_h (FMA)
 *   result = H + (R*p + (Tl[j] + (c1h*R)_l))
 * Lanes with x = +-0, +-inf or x < 0 (vfpclasspd 0x5e; NaN and denormals stay in the vector path) are
 * recomputed by a scalar routine, which also returns the error code and raises the FP flags.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_log28_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <math.h>
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../log2_d.h"

#if defined(__FAST_MATH__)
#error "svml_log2_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_log2_d_ha_constants.inc"

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

/* clears the low 32 bits (the assembly's "movl $0, -8(%rsp)" on a stored double) */
static inline double
hi32(double d)
{
    return bitsd(dbits(d) & 0xffffffff00000000ULL);
}

/*
 * __svml_dlog2_ha_cout_rare_internal: a libm-style scalar log2. Four paths by the high word of x:
 * x = 1 exactly, three polynomial paths for x within 2^-21, 2^-13 and 2^-6 of 1 (Q3, Q2, Q1, with
 * 32-bit splits for extra precision), and the general path with a 256-entry table of single-precision
 * reciprocals of the top 8 mantissa bits and {low, high} -log2(rcp). Statement order mirrors the assembly.
 */
int
svml_dlog2_ha_cout_rare(const double *px, double *pr)
{
    volatile double vzero = 0.0;  /* x - 0.0: quiets signaling NaNs (and raises invalid) like the assembly */
    double xd = *px - vzero;
    uint32_t hi = (uint32_t)(dbits(xd) >> 32);
    const uint32_t lo = (uint32_t)dbits(xd);

    if (hi >= 0x7ff00000u) {  /* negative, inf or NaN */
        const uint32_t ahi = hi & 0x7fffffffu;
        if (ahi > 0x7ff00000u || (ahi == 0x7ff00000u && lo != 0)) {  /* NaN */
            *pr = xd * sOnes[0];
            return 0;
        }
        if (!(hi & 0x80000000u)) {  /* +inf */
            *pr = xd;
            return 0;
        }
        if ((ahi | lo) != 0) {  /* x < 0 or -inf: NaN, invalid */
            volatile double inf = sInfs[0];
            *pr = copysign(inf * sZeros[0], -1.0);  /* x86's default NaN (negative) on every architecture */
            return 1;
        }
        volatile double z = sZeros[0];  /* -0: -inf, divide-by-zero */
        *pr = sOnes[1] / z;
        return 2;
    }

    if (hi - 0x3fef8000u < 0xd000u) {  /* x in [0.984375, 1.01953125): near 1 */
        if (hi - 0x3feffff0u >= 24u) {  /* Q1: degree-11 polynomial with 32-bit splits */
            const double t = xd + sOnes[1];
            const double t2 = t * t;
            double a = sQ1[9] * t2;
            double b = sQ1[8] * t2;
            const double t32p = sTwo32P * t;
            const double t32 = t * sTwo32;
            a = a + sQ1[7];
            b = b + sQ1[6];
            a = a * t2;
            b = b * t2;
            a = a + sQ1[5];
            b = b + sQ1[4];
            a = a * t2;
            b = b * t2;
            a = a + sQ1[3];
            b = b + sQ1[2];
            a = a * t2;
            b = b * t;
            a = a + b;                          /* A */
            const double a32p = sTwo32P * a;
            const double a32 = sTwo32 * a;
            const double th = t32p - t32;       /* t rounded to 21 bits */
            const double tl = t - th;
            const double ah = a32p - a32;
            double al = a - ah;
            double bb = sQ1[11] + ah;
            al = al + sQ1[1];
            double c = bb * tl;
            al = al * t;
            bb = bb * th;
            double ts = th + tl;
            c = c + al;
            double s = c + bb;
            const double s32p = sTwo32P * s;
            s = s * sTwo32;
            const double sh = s32p - s;
            c = c - sh;
            double d = sQ1[10] + sh;
            c = c + bb;
            c = c + sQ1[0];
            double f = d + c;
            const double f32p = sTwo32P * f;
            f = f * sTwo32;
            const double fh = f32p - f;
            ts = ts * fh;
            d = d - fh;
            c = c + d;
            c = c * t;
            *pr = ts + c;
            return 0;
        }
        const uint32_t top = (hi << 11) | (lo >> 21);
        if (top + 0x80000002u >= 3u) {  /* Q2 */
            const double t = xd + sOnes[1];
            const double t2 = t * t;
            double a = sQ2[2] * t2;
            double b = sQ2[3] * t2;
            a = a + sQ2[0];
            b = b + sQ2[1];
            a = a * t;
            b = b * t2;
            b = b + a;
            const double th = hi32(t);
            double s = t - th;
            s = th + s;
            s = s * sQ2[4];
            *pr = s + b;
            return 0;
        }
        if (((hi - 0x3ff00000u) | lo) == 0) {  /* x = 1 */
            *pr = sZeros[0];
            return 0;
        }
        const double t = xd + sOnes[1];  /* Q3 */
        double a = sQ3[2] * t;
        const double b = sQ3[3] * t;
        a = a + sQ3[1];
        a = a * t;
        a = a + sQ3[0];
        a = a * t;
        *pr = a + b;
        return 0;
    }

    /* general path */
    int e = (int)(hi >> 20) - 1023;
    uint32_t mhi = hi & 0xfffffu;
    if (hi < 0x00100000u) {  /* zero or denormal */
        if ((hi | lo) == 0) {  /* +0: -inf, divide-by-zero */
            volatile double z = sZeros[0];
            *pr = sOnes[1] / z;
            return 2;
        }
        xd = xd * sTwo55;
        mhi = (uint32_t)(dbits(xd) >> 32);
        e = (int)(mhi >> 20) - 1078;
        hi = mhi;
    }
    const uint32_t i = (hi >> 12) & 0xff;  /* top 8 mantissa bits */
    const double y = bitsd(((uint64_t)(mhi | 0x3ff00000u) << 32) | (uint32_t)dbits(xd));  /* in [1, 2) */
    const double rcp = (double)sRcp256[i];
    const double ed = (double)e;
    const double yh = hi32(y);
    double v = y - yh;
    v = v * rcp;
    double u = rcp * yh;
    u = u + sOnes[1];  /* rcp*yh - 1, exact */
    const double r = u + v;
    double up = u * sP[6];
    v = v * sP[6];
    const double up32 = sTwo32 + up;
    const double r2 = r * r;
    double a = sP[5] * r2;
    double b = sP[4] * r2;
    a = a + sP[3];
    b = b + sP[2];
    a = a * r2;
    b = b * r2;
    a = a + sP[1];
    b = b + sP[0];
    a = a * r2;
    b = b * r;
    a = a + b;
    const double eh = ed + sLog2Tab[2 * i + 1];
    v = v + a;
    const double uph = up32 - sTwo32;  /* u*P6 rounded to 20 fractional bits */
    up = up - uph;
    const double h = uph + eh;
    up = up + sLog2Tab[2 * i];
    v = v + up;
    *pr = v + h;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC

__attribute__((target("avx512f,avx512dq"))) __m512d
svml_log28_ha(__m512d x)
{
    const __m512d m = _mm512_getmant_round_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_nan, SAE);  /* [1,2), NaN for x < 0 */
    __m512d e = _mm512_getexp_round_pd(x, SAE);
    const __m512d rcp0 = _mm512_rcp14_pd(m);
    const __mmask8 special = _mm512_fpclass_pd_mask(x, 0x5e);  /* +-0, +-inf, negative */
    const __m512d rcp = _mm512_roundscale_round_pd(rcp0, 0x58, SAE);  /* RN to 5 fractional bits */
    const __m512d r = _mm512_fmsub_round_pd(rcp, m, V1(vOne), RN);   /* R = rcp*m - 1 */
    const __mmask8 lt34 = _mm512_cmp_round_pd_mask(rcp, V1(vThree4), _CMP_LT_OQ, SAE);
    const __m512d p910 = _mm512_fmadd_round_pd(V1(vC10), r, V1(vC9), RN);
    const __m512d p78 = _mm512_fmadd_round_pd(V1(vC8), r, V1(vC7), RN);
    const __m512d p56 = _mm512_fmadd_round_pd(V1(vC6), r, V1(vC5), RN);
    const __m512d r2 = _mm512_mul_round_pd(r, r, RN);
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(rcp), 48);
    const __m512d p710 = _mm512_fmadd_round_pd(r2, p910, p78, RN);
    e = _mm512_mask_add_round_pd(e, lt34, e, V1(vOne), RN);          /* E + 1 if rcp < 0.75 */
    const __m512d th = _mm512_permutex2var_pd(_mm512_loadu_pd(vTh), idx, _mm512_loadu_pd(vTh + 8));
    const __m512d tl = _mm512_permutex2var_pd(_mm512_loadu_pd(vTl), idx, _mm512_loadu_pd(vTl + 8));
    const __m512d r4 = _mm512_mul_round_pd(r2, r2, RN);
    const __m512d p34 = _mm512_fmadd_round_pd(V1(vC4), r, V1(vC3), RN);
    const __m512d kt = _mm512_add_round_pd(e, th, RN);               /* K + Th */
    const __m512d p36 = _mm512_fmadd_round_pd(r2, p56, p34, RN);
    const __m512d p12 = _mm512_fmadd_round_pd(V1(vC2), r, V1(vC1L), RN);
    __m512d p = _mm512_fmadd_round_pd(r4, p710, p36, RN);
    const __m512d h = _mm512_fmadd_round_pd(r, V1(vC1H), kt, RN);    /* K + Th + (c1h*R)_h */
    p = _mm512_fmadd_round_pd(r2, p, p12, RN);                        /* polynomial */
    __m512d crh = _mm512_sub_round_pd(h, kt, RN);                     /* (c1h*R)_h */
    crh = _mm512_mask_xor_pd(crh, special, crh, crh);
    const __m512d crl = _mm512_fmsub_round_pd(r, V1(vC1H), crh, RN);  /* (c1h*R)_l */
    const __m512d lo = _mm512_add_round_pd(tl, crl, RN);
    p = _mm512_fmadd_round_pd(r, p, lo, RN);
    __m512d res = _mm512_add_round_pd(p, h, RN);

    if (special) {
        double xin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (special & (1u << i)) {
                svml_dlog2_ha_cout_rare(&xin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

__attribute__((target("avx512f,avx512dq"))) void
svml_log2_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_log28_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_log28_ha(_mm512_mask_loadu_pd(_mm512_set1_pd(1.0), m, x + i)));
    }
}

#endif /* __x86_64__ */
