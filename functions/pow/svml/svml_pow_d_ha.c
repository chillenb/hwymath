/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision pow, high-accuracy variant
 * (svml_z0_pow_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), pow(x, y) = 2^(y*log2(x)):
 *   log2(x): m = getmant(x) in [1/2, 1), E = getexp(x) (exact, also for denormals)
 *     rcp = rcp14(m) rounded to a multiple of 1/32 (in (1, 2]); its top 5 mantissa bits index 32-entry
 *     tables Th + Tl = -log2(rcp) (+1 when rcp >= 1.5, which folds the binade of rcp into the index)
 *     R = (rcp*m - 1)/2 (exact, |R| < 2^-6), log2(x) = E + (rcp < 1.5) + Th + Tl + log2(1 + 2R)
 *     log2(1 + 2R) = C1H*R + R*(A0 + A1 R + ... + A9 R^9), C1H = 2/ln2 rounded, A0, A1 tiny corrections;
 *     the leading part E + Th + C1H*R is evaluated with error-free steps, giving High + Tll ~ log2(x)
 *   z = y*log2(x): P = RZ(y*High), Zl = RZ(y*High - P), T = RZ(y*Tll + Zl), S = RZ(P + T),
 *     low = T - (S - P), so S + low ~ z (round toward zero on the leading products and sums)
 *   2^z: S = N + j/16 + r with j, N from RD(S + 1.5*2^48 + 1023*16) (exponent bias in the low bits) and
 *     r = vreducepd(S) = S - floor(16 S)/16 (rounded down), Re = r + low
 *     2^z = 2^N * The[j] * (1 + Re*(B1 + B2 Re + ... + B6 Re^5) + Tle[j]/The[j])
 * Callout lanes: x <= 0 (or -0), x inf/NaN, y inf/NaN, and |S| > 1021.5 or NaN (overflow, underflow,
 * denormal results). They are recomputed by a scalar routine with its own algorithm (3-step reciprocal
 * table reduction for log2, 2^(j/128) table for exp2) and special-case logic.
 *
 * Build: any x86-64 C compiler with AVX-512F/DQ intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_pow8_ha() needs AVX-512 at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../pow_d.h"

#if defined(__FAST_MATH__)
#error "svml_pow_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_pow_d_ha_constants.inc"

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

/* the assembly reads and writes parts of doubles through the stack: high 32 bits, top 16 bits */
static inline uint32_t hi32(double d) { return (uint32_t)(dbits(d) >> 32); }
static inline uint32_t lo32(double d) { return (uint32_t)dbits(d); }
static inline uint32_t top16(double d) { return (uint32_t)(dbits(d) >> 48); }
static inline double
with_top16(double d, uint32_t t)
{
    return bitsd((dbits(d) & 0x0000ffffffffffffULL) | ((uint64_t)(t & 0xffff) << 48));
}

/* x86 32-bit shift (count taken modulo 32) */
static inline uint32_t shl32(uint32_t v, uint32_t c) { return v << (c & 31); }

/* Dekker split with the 2^27+1 splitter: hi has at most 26 significant bits, hi + lo = a */
static inline double
split_hi(double a)
{
    const double s = sSplitter * a;
    const double t = s - a;
    return s - t;
}

/*
 * __svml_dpow_ha_cout_rare_internal: scalar pow for the lanes flagged by the vector path.
 * Statement order mirrors the assembly (labels in comments); values that only exist to raise floating-
 * point exception flags are computed through volatile variables so that they happen at run time.
 */
int
svml_dpow_ha_cout_rare(const double *px, const double *py, double *pr)
{
    volatile double one = sOne;  /* x*1, y*1: quiet signaling NaNs (raising invalid) like the assembly */
    const double x1 = *px * one;
    const double y1 = *py * one;
    int ret = 0;

    const uint32_t ex = (top16(x1) & 0x7ff0) >> 4;  /* r13d */
    const uint32_t ey = (top16(y1) & 0x7ff0) >> 4;  /* edx */
    const uint32_t sx = top16(x1) >> 15;            /* r12d */
    const uint32_t sy = top16(y1) >> 15;            /* bl */
    const int xmant0 = (hi32(x1) & 0xfffff) == 0 && lo32(x1) == 0;  /* r11b */
    const uint32_t yhi20 = hi32(y1) & 0xfffff;
    const uint32_t ylo = lo32(y1);
    int ymant0 = yhi20 == 0 && ylo == 0;  /* r10d */
    const uint32_t ywords = (ey << 20) | yhi20;
    const int y_finite = ey != 2047;  /* esi */
    uint32_t yint;  /* r9d: 0 = y not an integer, 1 = odd integer, 2 = even integer (or 0) */
    int x_finite;   /* cl */
    double res;     /* xmm2 */

    /* classify y */
    if ((ylo | ywords) == 0) { /* .LBL_2_24: y = +-0 */
        yint = 2;
        goto L25;
    }
    if (ey < 1023 || !y_finite) { /* .LBL_2_23: |y| < 1, inf, NaN */
        yint = 0;
        goto L25;
    }
    {
        const uint32_t e = ey - 1023;
        if (e <= 20) {
            if ((shl32(shl32(ywords, e), 12) | ylo) != 0) {
                yint = 0;
            }
            else { /* .LBL_2_15: the units bit decides odd/even */
                yint = (~shl32(ywords, ey - 1012) >> 31) + 1;
            }
        }
        else if (e < 53) {
            const uint32_t t = shl32(ylo, ey - 1012);
            if (t & 0x7fffffff) {
                yint = 0;
            }
            else {
                yint = (~t >> 31) + 1;
            }
        }
        else { /* .LBL_2_20 */
            yint = 2;
        }
    }
    /* .LBL_2_21 */
    if (sx == 0 && ex == 1023) { /* .LBL_2_57: x in [1, 2) */
        if (xmant0) {
            goto L55;
        }
        goto L30;
    }
    goto L29;

L25:
    if (sx == 0 && ex == 1023) { /* .LBL_2_74 */
        if (xmant0) {
            goto L55;
        }
        if (ey != 0) {
            goto L30;
        }
        if (!ymant0) {
            goto L30;
        }
        goto L55;
    }
    /* .LBL_2_27 */
    if (ey == 0 && ymant0) {
        goto L55;
    }

L29:
    if (ex == 2047) {
        x_finite = 0;
        goto L32;
    }
L30:
    x_finite = 1;
L32:
    if (!(xmant0 | x_finite)) { /* x NaN */
        goto L54;
    }
    ymant0 |= y_finite;
    if (!ymant0) { /* y NaN */
        goto L54;
    }
    res = sZero;
    if (x1 == sZero) { /* .LBL_2_71 */
        const double xx = x1 * x1;
        if (sy == 0) { /* .LBL_2_73: +-0 */
            res = (yint & sx ? sMOne : sOne) * xx;
        }
        else { /* pole: +-inf, raises divide-by-zero */
            res = (yint & sx ? sMOne : sOne) / xx;
            ret = 1;
        }
        goto L56;
    }
    if (x1 == sMOne) { /* .LBL_2_68 */
        if (yint != 0 || !y_finite) { /* .LBL_2_70 */
            res = (yint & 1) ? sMOne : sOne;
            goto L56;
        }
    }
    /* .LBL_2_36 */
    if (!x_finite || !y_finite) { /* .LBL_2_47 */
        if (ex < 1023) {
            if (sy) {
                res = y1 * y1;
            }
            goto L56;  /* else res = +0 */
        }
        /* .LBL_2_51 */
        if (sy) {
            res = res * ((yint & sx) ? sMOne : sOne);
        }
        else { /* .LBL_2_53 */
            double t = x1 * x1;
            t = t * y1;
            res = t * ((yint & sx) ? sMOne : sOne);
        }
        goto L56;
    }
    if (!(x1 > sZero) && yint == 0) { /* .LBL_2_46: x < 0, y not an integer: NaN, raises invalid */
        volatile double z = res;
        res = z / z;
        ret = 1;
        goto L56;
    }

    { /* .LBL_2_40: the general case */
        const double sign = (yint & sx) ? sMOne : sOne;  /* xmm0 */
        double ax = bitsd(dbits(x1) & 0x7fffffffffffffffULL);
        int32_t escale = 0;
        if (ex == 0) { /* denormal x */
            ax = ax * sTwo200;
            escale = -200;
        }
        /* .LBL_2_43: |x| = 2^k * X1, X1 in [1, 2); k is one larger when X1 >= 1.484375 (Rcp1 index >= 16) */
        int32_t k = (int32_t)(hi32(ax) - 0x3fe7c000u) >> 20;
        k += escale;
        const double m = with_top16(ax, (top16(ax) & 0x800f) + 0x3ff0);  /* X1 */
        const uint32_t i1 = ((hi32(m) & 0xfc000) + 0x4000) >> 15;
        double rcp = sRcp1[i1];  /* xmm5 */
        double xr = rcp * m;     /* xmm6: X2 = Rcp1*X1 */
        double sm = sSplitter * m;
        double hh = (double)k;   /* xmm7: k + sum of log2(1/Rcp) (high parts) */
        hh = hh + sLog1[2 * i1];
        const uint32_t i2 = ((hi32(xr) & 0xfc00) + 0x400) >> 11;
        const double smm = sm - m;
        const double rcp2 = sRcp2[i2];
        xr = xr * rcp2;
        hh = hh + sLog2[2 * i2];
        rcp = rcp * rcp2;
        const double mh = sm - smm;  /* split of X1 */
        const uint32_t i3 = ((hi32(xr) & 0xff0) + 16) >> 5;
        const double rcp3 = sRcp3[i3];
        rcp = rcp * rcp3;  /* RcpC = Rcp1*Rcp2*Rcp3C, exact */
        const double ml = m - mh;
        hh = hh + sLog3[2 * i3];
        xr = xr * rcp3;   /* X1*RcpC, rounded */
        double cq = xr - sCHi;
        const double rh = split_hi(rcp);
        const double lo1 = sLog1[2 * i1 + 1];
        const double rl = rcp - rh;
        const double lo2 = sLog2[2 * i2 + 1];
        const double lo3 = sLog3[2 * i3 + 1];
        double err = mh * rh;  /* xmm2: X1*RcpC - xr, exactly (Dekker) */
        const double t5 = ml * rh;
        err = err - xr;
        double h = hh;
        err = err + t5;
        h = h + cq;           /* HH */
        err = err + mh * rl;
        err = err + ml * rl;
        double hl = hh - h;
        hl = hl + cq;         /* HL */
        const double HH = h;

        const uint32_t ez = ((top16(HH) & 0x7ff0) >> 4) + ey;
        if (ez >= 2057) { /* .LBL_2_67: |y*log2|x|| huge: overflow or underflow */
            volatile double v = ((sy ^ (top16(HH) >> 15)) & 1) ? sTiny : sHuge;  /* raises the flags */
            res = v * v;
            res = res * sign;
            goto L56;
        }
        if (ez <= 1984) { /* |y*log2|x|| tiny: +-(1 + tiny), raises inexact */
            volatile double t = sOne;
            res = t + sTiny;
            res = res * sign;
            goto L56;
        }

        /* .LBL_2_58: log2|x| = HH + HL + low parts, polynomial in cq */
        cq = cq + err;
        double lo = lo2 + lo3;
        err = err + hl;
        double p = sL4 * cq;
        lo = lo1 + lo;
        p = p + sL3;
        lo = lo + err;
        p = p * cq;
        p = p + sL2;
        p = p * cq;
        p = p + sL1;
        p = p * cq;
        double H1 = HH + p;
        double e1 = HH - H1;
        e1 = e1 + p;
        double H = lo + H1;  /* log2|x| = H + (e2 + e1) */
        double e2 = H1 - H;
        e2 = e2 + lo;
        const double Hs = sSplitter * H;
        double L = e2 + e1;
        const double y = y1;
        const double Hsh = Hs - H;
        L = L * y;
        const double ys = y * sSplitter;
        const double Hh = Hs - Hsh;
        const double Hl = H - Hh;
        const double ysy = ys - y;
        const double yh = ys - ysy;
        const double yl = y - yh;
        const double pp = yh * Hh;  /* exact */
        double z = yh * Hl;
        double sh = sShifter + pp;  /* integer bits: 128*N + j */
        z = z + Hl * yl;
        z = z + Hh * yl;
        z = z + L;
        const double nj = sh - sShifter;
        const double zr = pp - nj;
        const uint32_t ib = lo32(sh);
        const uint32_t j = ib & 127;
        double r = zr + z;  /* y*log2|x| = N + j/128 + r */
        double q = sE5 * r;
        const double thi = sExp[2 * j];
        const int32_t n = (int32_t)ib >> 7;
        q = q + sE4;
        q = q * r;
        q = q + sE3;
        q = q * r;
        q = q + sE2;
        q = q * r;
        q = q + sE1;
        q = q * r;
        q = q * thi;
        q = q + sExp[2 * j + 1];
        const double s = q + thi;  /* 2^(j/128 + r) */
        const uint32_t ts = top16(s);
        const int32_t e = n + (int32_t)((ts & 0x7ff0) >> 4) - 1023;
        if (e >= 1024) { /* .LBL_2_66: overflow */
            volatile double v = sHuge;
            res = v * v;
            res = sign * res;
            goto L56;
        }
        if (e >= -1022) { /* normal result: set the exponent */
            res = with_top16(s, (ts & 0x800f) | ((uint32_t)(e + 1023) & 2047) << 4);
            res = res * sign;
            goto L56;
        }
        if (e >= -1032) { /* .LBL_2_61: denormal result, one rounding */
            const uint32_t be = ((uint32_t)n + 1223) & 2047;  /* 2^(N+200) */
            double sl = thi - s;
            sl = q + sl;
            const double s_h = split_hi(s);
            const double s_l = s - s_h;
            const double scale = bitsd((uint64_t)((0x3ff0 & 0x800f) | be << 4) << 48);
            const double a = (s_l + sl) * scale;
            const double b = s_h * scale;
            double t = sign * sTwoM200;
            const double c = a + b;
            volatile double tiny = sTiny;  /* raises underflow at run time */
            const double tt = tiny * tiny;
            t = t * c;
            res = t + tt;
            goto L56;
        }
        if (e >= -1084) { /* .LBL_2_63 */
            const uint32_t be = ((uint32_t)n + 1223) & 2047;
            const double scale = bitsd((uint64_t)((0x3ff0 & 0x800f) | be << 4) << 48);
            double t = sign * sTwoM200;
            const double a = s * scale;
            t = t * a;
            volatile double tiny = sTiny;
            const double tt = tiny * tiny;
            res = t - tt;
            goto L56;
        }
        { /* .LBL_2_65: underflow to +-0 */
            volatile double tiny = sTiny;
            const double tt = tiny * tiny;
            res = sign * tt;
            goto L56;
        }
    }

L54: /* NaN operand */
    res = x1 + y1;
    goto L56;

L55: /* pow(1, y) = pow(x, +-0) = 1 */
    {
        volatile double s = x1 + y1;  /* may raise inexact, like the assembly */
        const double v = bitsd((dbits(s) & 0x8000000000000000ULL) | dbits(sOne));
        res = v * v;
    }

L56:
    *pr = res;
    return ret;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define I1(c) _mm512_set1_epi64((long long)(c))
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define RZ (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)
#define RD (_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC

#define AVX512 __attribute__((target("avx512f,avx512dq")))

/* 32-entry table lookup: vpermt2pd on the low 4 index bits, blend on bit 4 */
static AVX512 inline __m512d
lookup32(const double *t, __m512i idx, __mmask8 hi)
{
    const __m512d a = _mm512_permutex2var_pd(_mm512_loadu_pd(t), idx, _mm512_loadu_pd(t + 8));
    const __m512d b = _mm512_permutex2var_pd(_mm512_loadu_pd(t + 16), idx, _mm512_loadu_pd(t + 24));
    return _mm512_mask_blend_pd(hi, a, b);
}

AVX512 __m512d
svml_pow8_ha(__m512d x, __m512d y)
{
    /* log2(x) */
    const __m512d mant = _mm512_getmant_round_pd(x, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_nan, SAE);  /* [1/2, 1) */
    __m512d expon = _mm512_getexp_round_pd(x, SAE);
    const __m512d rcp0 = _mm512_rcp14_pd(mant);
    const __mmask8 x_special = _mm512_fpclass_pd_mask(x, 0xdf);  /* x <= 0 (not +denormal), inf, NaN */
    const __m512d rcp = _mm512_roundscale_round_pd(rcp0, 0x58, SAE);  /* to 1/32, round to nearest */
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(rcp), 47);
    const __mmask8 lt15 = _mm512_cmp_round_pd_mask(rcp, V1(vRcpLimit), _CMP_LT_OQ, SAE);
    const __mmask8 hi16 = _mm512_movepi64_mask(_mm512_slli_epi64(idx, 59));  /* index bit 4 */
    const __m512d th = lookup32(vLogTh, idx, hi16);
    expon = _mm512_mask_add_round_pd(expon, lt15, expon, V1(vOne), RN);  /* +1 if rcp < 1.5 */
    const __m512d tl = lookup32(vLogTl, idx, hi16);
    const __m512d mh = _mm512_mul_round_pd(mant, V1(vHalf), RN);
    const __mmask8 y_special = _mm512_fpclass_pd_mask(y, 0x99);  /* y inf or NaN */
    const __m512d R = _mm512_fmsub_round_pd(rcp, mh, V1(vHalf), RN);  /* (rcp*m - 1)/2, exact */

    const __m512d R2 = _mm512_mul_round_pd(R, R, RN);
    __m512d p910 = _mm512_fmadd_round_pd(R, V1(vA9), V1(vA8), RN);
    const __m512d p78 = _mm512_fmadd_round_pd(R, V1(vA7), V1(vA6), RN);
    const __m512d eth = _mm512_add_round_pd(th, expon, RN);  /* Expon + Th */
    const __m512d p710 = _mm512_fmadd_round_pd(p910, R2, p78, RN);
    const __m512d p12 = _mm512_fmadd_round_pd(R, V1(vA1), V1(vA0), RN);
    const __m512d p56 = _mm512_fmadd_round_pd(R, V1(vA5), V1(vA4), RN);
    const __m512d p34 = _mm512_fmadd_round_pd(R, V1(vA3), V1(vA2), RN);
    const __m512d R4 = _mm512_mul_round_pd(R2, R2, RN);
    const __m512d p36 = _mm512_fmadd_round_pd(p56, R2, p34, RN);
    const __m512d p310 = _mm512_fmadd_round_pd(p710, R4, p36, RN);
    const __m512d p110 = _mm512_fmadd_round_pd(p310, R2, p12, RN);
    const __m512d tlp = _mm512_fmadd_round_pd(p110, R, tl, RN);  /* Tl + R*P1_10 */

    const __m512d s1 = _mm512_fmadd_round_pd(R, V1(vC1H), eth, RN);  /* Expon + Th + C1H*R */
    const __m512d rh = _mm512_sub_round_pd(s1, eth, RN);             /* (C1H*R)_high */
    const __m512d s2 = _mm512_fnmadd_round_pd(R, rh, s1, RN);        /* - (R*Rh)_high */
    const __m512d rl = _mm512_fmsub_round_pd(R, V1(vC1H), rh, RN);   /* (C1H*R)_low */
    const __m512d rrh = _mm512_sub_round_pd(s1, s2, RN);             /* (R*Rh)_high */
    const __m512d rl2 = _mm512_fnmadd_round_pd(rl, R, rl, RN);       /* Rl - R*Rl */
    const __m512d rrl = _mm512_fmsub_round_pd(rh, R, rrh, RN);       /* (R*Rh)_low */
    const __m512d low1 = _mm512_sub_round_pd(rl2, rrl, RN);
    const __m512d low = _mm512_add_round_pd(tlp, low1, RN);
    const __m512d high = _mm512_add_round_pd(s2, low, RN);          /* log2(x), high part */

    /* y*log2(x), rounding toward zero */
    const __m512d P = _mm512_mul_round_pd(y, high, RZ);
    const __m512d tlh = _mm512_sub_round_pd(high, s2, RN);
    const __m512d zl = _mm512_fmsub_round_pd(y, high, P, RZ);       /* (y*High)_low */
    const __m512d tll = _mm512_sub_round_pd(low, tlh, RN);
    const __m512d T = _mm512_fmadd_round_pd(y, tll, zl, RZ);        /* y*Tll + Zl */
    const __m512d S = _mm512_add_round_pd(P, T, RZ);
    const __m512d sp = _mm512_sub_round_pd(S, P, RN);
    const __m512d sh = _mm512_add_round_pd(S, V1(vShifter), RD);

    /* 2^(S + low) */
    const __m512d red = _mm512_reduce_round_pd(S, 0x41, SAE);       /* S - floor(16 S)/16 */
    const __m512d lw = _mm512_sub_round_pd(T, sp, RN);
    const __m512d aS = _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(S), I1(vAbsMaskBits)));
    __m512d re = _mm512_add_round_pd(red, lw, RN);
    re = _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(re), I1(vRMaskBits)));  /* |Re| < 2 */
    const __m512i jdx = _mm512_castpd_si512(sh);
    const __m512d tle = _mm512_permutex2var_pd(_mm512_loadu_pd(vExpTle), jdx, _mm512_loadu_pd(vExpTle + 8));
    const __m512d the = _mm512_permutex2var_pd(_mm512_loadu_pd(vExpThe), jdx, _mm512_loadu_pd(vExpThe + 8));
    const __m512d re2 = _mm512_mul_round_pd(re, re, RN);
    const __mmask8 range = _mm512_cmp_round_pd_mask(aS, V1(vBound), _CMP_NLE_UQ, SAE);
    __m512d q = _mm512_fmadd_round_pd(re, V1(vB6), V1(vB5), RN);
    const __m512d q34 = _mm512_fmadd_round_pd(re, V1(vB4), V1(vB3), RN);
    const __m512d q12 = _mm512_fmadd_round_pd(re, V1(vB2), V1(vB1), RN);
    __m512i sc = _mm512_slli_epi64(jdx, 48);
    q = _mm512_fmadd_round_pd(q, re2, q34, RN);
    sc = _mm512_and_si512(sc, I1(vExpMaskBits));  /* 2^N */
    q = _mm512_fmadd_round_pd(q, re2, q12, RN);
    q = _mm512_fmadd_round_pd(q, re, tle, RN);     /* R*poly + Tle/The */
    q = _mm512_fmadd_round_pd(q, the, the, RN);    /* The + The*(...) */
    __m512d res = _mm512_mul_round_pd(q, _mm512_castsi512_pd(sc), RN);

    const unsigned callout = (unsigned)(x_special | y_special | range);
    if (callout) {
        double xin[8], yin[8], out[8];
        _mm512_storeu_pd(xin, x);
        _mm512_storeu_pd(yin, y);
        _mm512_storeu_pd(out, res);
        for (int i = 0; i < 8; i++) {
            if (callout & (1u << i)) {
                svml_dpow_ha_cout_rare(&xin[i], &yin[i], &out[i]);
            }
        }
        res = _mm512_loadu_pd(out);
    }
    return res;
}

AVX512 void
svml_pow_ha_array(const double *x, const double *y, double *z, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(z + i, svml_pow8_ha(_mm512_loadu_pd(x + i), _mm512_loadu_pd(y + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(z + i, m, svml_pow8_ha(_mm512_maskz_loadu_pd(m, x + i), _mm512_maskz_loadu_pd(m, y + i)));
    }
}

#endif /* __x86_64__ */
