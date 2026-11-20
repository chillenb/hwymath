/* Scalar reference of VRCP14PD/VRSQRT14PD (double precision), after Intel's RECIP14.c (Copyright (c) 2015
 * Intel Corporation, BSD-3-Clause), for normal inputs whose results are normal, +-0, +inf and NaN.
 * Copied from functions/asinh/tests/approx14_ref.h (branch port/invhyperbolic).
 * test_approx14 checks it against the hardware instructions (or Intel SDE) and the Highway emulation. */
#ifndef HWYINVTRIG_APPROX14_REF_H
#define HWYINVTRIG_APPROX14_REF_H

#include <math.h>
#include <stdint.h>

#include "functions/asin/hwy/approx14_tables.h"
#include "hwymath/test_util.h"

static inline double
ref_rcp14(double x)
{
    const uint64_t u = dbits(x), s = u & 0x8000000000000000ULL, e = (u >> 52) & 0x7ff, m = u & 0x000fffffffffffffULL;
    if (isnan(x)) return bitsd(u | 0x0008000000000000ULL);
    if (isinf(x)) return bitsd(s);
    if (e == 0 && m == 0) return bitsd(s | 0x7ff0000000000000ULL);
    if (m == 0) return bitsd(s | ((0x7feULL - e) << 52));  /* powers of two (normal results) */
    const double xm = bitsd(m | 0x3ff0000000000000ULL);    /* in (1, 2) */
    const unsigned i = (unsigned)floor((xm - 1.0) * 64);
    const double y = 1.0 + 1.0 / 128 + i / 64.0;
    const double xt = bitsd(dbits(xm) & 0xfffffff000000000ULL);
    const double r = (double)kRcp14Coeff[2 * i + 1] - 256.0 * kRcp14Coeff[2 * i] * (xt - y);
    const uint64_t r12 = (dbits(r) & 0xfffffff000000000ULL) - 0x0120000000000000ULL;  /* in (0.5, 1) */
    return bitsd(s | (r12 & 0x000fffffffffffffULL) | ((0x7fdULL - e) << 52));
}

static inline double
ref_rsqrt14_14(double x)  /* x in [1, 4) */
{
    if (x == 1.0) return 1.0;
    unsigned i;
    double y, c, dd;
    if (x < 2.0) {
        i = (unsigned)floor((x - 1.0) * 32);
        y = 1.0 + 1.0 / 64 + i / 32.0;
        c = kRsqrt14Coeff[4 * i + 1];
        dd = 256.0 * kRsqrt14Coeff[4 * i];
    } else {
        i = (unsigned)floor((x - 2.0) * 16);
        y = 2.0 + 1.0 / 32 + i / 16.0;
        c = kRsqrt14Coeff[4 * i + 3];
        dd = 128.0 * kRsqrt14Coeff[4 * i + 2];
    }
    const double xt = bitsd(dbits(x) & 0xffffffe000000000ULL);
    const double r = c - dd * (xt - y);
    return bitsd((dbits(r) & 0xfffffff000000000ULL) - 0x0130000000000000ULL);
}

static inline double
ref_rsqrt14(double x)
{
    const uint64_t u = dbits(x);
    if (isnan(x)) return bitsd(u | 0x0008000000000000ULL);
    if (x == 0.0) return bitsd((u & 0x8000000000000000ULL) | 0x7ff0000000000000ULL);
    if (x < 0.0) return bitsd(0xfff8000000000000ULL);
    if (isinf(x)) return 0.0;
    const int64_t ex = (int64_t)(u >> 52) - 1023;  /* normal x */
    const uint64_t m = u & 0x000fffffffffffffULL;
    const int64_t n = ex >> 1;  /* floor(ex/2) */
    const double r = ref_rsqrt14_14(bitsd(m | ((ex & 1) ? 0x4000000000000000ULL : 0x3ff0000000000000ULL)));
    return bitsd(dbits(r) - ((uint64_t)n << 52));
}

#endif
