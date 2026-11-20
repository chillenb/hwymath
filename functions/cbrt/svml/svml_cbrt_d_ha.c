/*
 * C/intrinsics translation of Intel SVML's AVX-512 double-precision cbrt, high-accuracy variant
 * (svml_z0_cbrt_d_ha.s from https://github.com/numpy/SVML), keeping its algorithm, constants,
 * operation order and rounding modes, so that results should be bit-identical.
 *
 * Original: Copyright (C) 2023 Intel Corporation, SPDX-License-Identifier: BSD-3-Clause.
 * This translation is distributed under the same license.
 *
 * Algorithm (vector path, 8 lanes), x = 2^E * m with m = getmant(x) in [1, 2) (sign of x), E = getexp(x):
 *   rcp  = 1/m to 4 fractional bits: vrcp14pd (~14-bit reciprocal) rounded to a multiple of 1/16, k/16
 *          with k = 8..16 (so rcp is exactly representable and the table index is its exponent/top bits)
 *   R    = m*rcp - 1                     (one FMA; |R| <~ 1/32)
 *   K    = floor(E/3), j = E - 3K         (as floor((E + 1.5*2^52)/3*(1+eps) - 2^51), exact for all E)
 *   Sh + Sl = 2^(j/3), Th + Tl = cbrt(1/rcp) = cbrt(16/k)  (8-entry vpermpd / 16-entry vpermt2pd tables)
 *   cbrt(x) = 2^K * Sh*Th * cbrt(1 + R)
 *           = 2^K * [ ShTh + (ShTh_lo + Sh*Tl + Sl*Th) + ShTh*R*P(R) ],  P(R) = (cbrt(1+R) - 1)/R
 *   with ShTh = RN(Sh*Th) and its rounding error ShTh_lo from an FMA, and the sign of x ORed in.
 * There is no scalar callout in the vector path: getexp/getmant handle denormals, and 0, inf and NaN
 * come out right by themselves (0 -> 2^-inf scaling, inf -> 2^inf scaling, NaN propagates).
 *
 * The scalar routine __svml_dcbrt_ha_cout_rare_internal (VML's scalar fallback; not called here) uses a
 * 32-entry table rc ~ m^(-1/3) and cbrt(m) = m*rc^2*(1 - d)^(-2/3) with d = 1 - m*rc^3, evaluated in
 * double-double arithmetic with Dekker splits.
 *
 * Build: any x86-64 C compiler with AVX-512F intrinsics (GCC, Clang, ICX). Do not use -ffast-math.
 * Only svml_cbrt8_ha() needs AVX-512F at run time; the scalar routine is plain C.
 */
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <immintrin.h>
#endif

#include "../cbrt_d.h"

#if defined(__FAST_MATH__)
#error "svml_cbrt_d_ha.c relies on exact IEEE arithmetic; do not compile with -ffast-math"
#endif
/* The scalar routine uses separate multiplies and adds; contracting them into FMAs changes results. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "svml_cbrt_d_ha_constants.inc"

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
 * __svml_dcbrt_ha_cout_rare_internal: scalar cbrt. x = 2^(3q + r - 1023) * mx, mx in [1, 2);
 * cbrt(x) = 2^(q - 341) * 2^(r/3) * cbrt(mx), cbrt(mx) = mx*rc^2 * (1 - d)^(-2/3), d = 1 - mx*rc^3.
 * Statement order (and the order of the operands of every addition) mirrors the assembly.
 */
int
svml_dcbrt_ha_cout_rare(const double *px, double *pr)
{
    const double x = *px;
    const uint64_t xb = dbits(x);
    const uint32_t ex = (uint32_t)(xb >> 52) & 0x7ff;

    if (ex == 0x7ff) {  /* inf, NaN */
        *pr = x + x;
        return 0;
    }
    if (x == sZero) {  /* +-0 */
        *pr = x * sOne;
        return 0;
    }

    const double sgn = sSign[xb >> 63];
    double ax = bitsd(xb & 0x7fffffffffffffffULL);
    uint32_t adj = 0;
    if (ex == 0) {  /* denormal: scale by 2^300, result by 2^-100 */
        ax = ax * sTwo300;
        adj = 100;
    }
    const uint32_t e = (uint32_t)(dbits(ax) >> 52) & 0x7ff;
    const uint32_t q = e / 3;
    const uint32_t r = e - 3 * q;
    const uint32_t be = (q + 682 - adj) & 2047;  /* biased exponent of 2^(q - 341 - adj) */
    const double mx = bitsd((dbits(ax) & 0x800fffffffffffffULL) | 0x3ff0000000000000ULL);

    /* mx = mh + ml, mh with 10 significant bits */
    double v = mx * sSplit43;
    double t = v - mx;
    const double mh = v - t;
    const double ml = mx - mh;
    const uint32_t i = (uint32_t)(dbits(mx) >> 47) & 31;  /* top 5 mantissa bits */
    const double rc = sRcp[i];
    const double rc2 = rc * rc;
    const double a = mh * rc2;  /* exact */
    const double d = sOne - rc * a;
    const double b = ml * rc2;
    v = sSplit27 * d;
    double c = rc * b;
    c = -c;
    t = v - d;
    const double dh = v - t;
    const double dl = d - dh;
    const double el = c + dl;  /* d = 1 - rc^3*mx ~ dh + el */
    const double z = dh + el;

    double p = sC12 * z;  /* p = z*(c2 + c3 z + ... + c12 z^10) */
    p = p + sC11;
    p = p * z;
    p = p + sC10;
    p = p * z;
    p = p + sC9;
    p = p * z;
    p = p + sC8;
    p = p * z;
    p = p + sC7;
    p = p * z;
    p = p + sC6;
    p = p * z;
    p = p + sC5;
    p = p * z;
    p = p + sC4;
    p = p * z;
    p = p + sC3;
    p = p * z;
    p = p + sC2;
    p = p * z;

    /* s + sl = 2/3 + p (two-sum, plus the low part of 2/3) */
    const double s = sTwoThirds + p;
    double u = sTwoThirds - s;
    double w = s + u;
    u = u + p;
    w = sTwoThirds - w;
    const double sl = u + w;
    v = sSplit27 * s;
    double sc = sl + sTwoThirdsLo;
    t = v - s;
    const double sh = v - t;
    const double sL = (s - sh) + sc;

    /* ph + pL = (dh + el)*(sh + sL) = d*(2/3 + ...) */
    const double P = dh * sh;
    double pl = sh * el;
    pl = pl + el * sL;
    pl = pl + dh * sL;
    v = sSplit27 * P;
    t = v - P;
    const double ph = v - t;
    const double pL = pl + (P - ph);

    /* cbrt(mx) ~ (mh + ml)*rc^2*(1 + ph + pL) = (a + b)*(1 + ph + pL) = B + Bl + A + ... */
    const double phb = ph * b;
    const double pha = ph * a;
    const double B = pha + a;
    double A = pL * b;
    A = A + phb;
    A = pL * a + A;
    const double aB = a - B;
    const double Bl = (pha + aB) + (a - (B + aB));
    v = sSplit27 * B;
    t = v - B;
    const double Bh = v - t;
    double L = (B - Bh) + A;
    L = L + b;
    L = Bl + L;

    /* * 2^(r/3) = hi2 + lo2, * 2^(q-341-adj), sign */
    const double hi2 = sCbrt2[2 * r];
    const double lo2 = sCbrt2[2 * r + 1];
    double lo = lo2 * L;
    lo = lo + lo2 * Bh;
    lo = lo + L * hi2;
    double res = Bh * hi2 + lo;
    res = res * bitsd((uint64_t)be << 52);
    *pr = sgn * res;
    return 0;
}

#if defined(__x86_64__)  /* the vector routines; the scalar callout above is portable C */

#define V1(c) _mm512_set1_pd(c)
#define RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)
#define SAE _MM_FROUND_NO_EXC

/* hides the table contents from the optimizer, which would otherwise replace the load of the all-zero
 * upper half of vTl by a zeroing instruction (same results, one extra instruction) */
static inline const double *
opaque(const double *p)
{
    __asm__("" : "+r"(p));
    return p;
}

__attribute__((target("avx512f"))) __m512d
svml_cbrt8_ha(__m512d x)
{
    const __m512d E = _mm512_getexp_round_pd(x, SAE);                                   /* floor(log2|x|) */
    const __m512d m = _mm512_getmant_round_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_src, SAE); /* +-[1, 2) */
    const __m512d rcp = _mm512_rcp14_pd(m);
    const __m512d ES = _mm512_add_round_pd(E, V1(vShifter), RN);                        /* 1.5*2^52 + E */
    const __m512i sgn = _mm512_and_si512(_mm512_castpd_si512(x), _mm512_set1_epi64((long long)vSignMaskBits));
    /* rcp to 4 fractional bits (imm 0x48: M = 4, round to nearest, suppress precision exception) */
    const __m512d rcpR = _mm512_roundscale_round_pd(rcp, 0x48, SAE);
    /* E/3 + 1/4 (+ tiny): floor gives floor(E/3) */
    const __m512d t = _mm512_fmsub_round_pd(V1(vOneThird), ES, V1(vTwo51), RN);
    const __m512d R = _mm512_fmsub_round_pd(m, rcpR, V1(vOne), RN);                     /* m*rcp - 1 */
    const __m512d K = _mm512_roundscale_round_pd(t, 0x09, SAE);                        /* floor */
    /* Th/Tl index: bits 49..52 of rcp = k - 8 (exponent LSB and top 3 mantissa bits) */
    const __m512i idx = _mm512_srli_epi64(_mm512_castpd_si512(rcpR), 49);
    /* 1.5*2^52 + E - 3K: j = E mod 3 in the low bits */
    const __m512d jS = _mm512_fnmadd_round_pd(V1(vThree), K, ES, RN);
    __m512d p = _mm512_fmadd_round_pd(V1(vC9), R, V1(vC8), RN);
    const __m512i jidx = _mm512_castpd_si512(jS);
    const __m512d Sh = _mm512_permutexvar_pd(jidx, _mm512_loadu_pd(vSh));
    const double *th = opaque(vTh), *tl = opaque(vTl);
    const __m512d Th = _mm512_permutex2var_pd(_mm512_loadu_pd(th), idx, _mm512_loadu_pd(th + 8));
    const __m512d Tl = _mm512_permutex2var_pd(_mm512_loadu_pd(tl), idx, _mm512_loadu_pd(tl + 8));
    const __m512d Sl = _mm512_permutexvar_pd(jidx, _mm512_loadu_pd(vSl));

    const __m512d ShTh = _mm512_mul_round_pd(Sh, Th, RN);                               /* (Sh*Th)_high */
    const __m512d R2 = _mm512_mul_round_pd(R, R, RN);
    const __m512d q45 = _mm512_fmadd_round_pd(V1(vC5), R, V1(vC4), RN);
    const __m512d q67 = _mm512_fmadd_round_pd(V1(vC7), R, V1(vC6), RN);
    p = _mm512_fmadd_round_pd(R2, p, q67, RN);
    __m512d lo = _mm512_fmsub_round_pd(Th, Sh, ShTh, RN);                               /* (Sh*Th)_low */
    p = _mm512_fmadd_round_pd(R2, p, q45, RN);
    lo = _mm512_fmadd_round_pd(Tl, Sh, lo, RN);                                         /* + Sh*Tl */
    lo = _mm512_fmadd_round_pd(Th, Sl, lo, RN);                                         /* + Sl*Th */
    const __m512d q23 = _mm512_fmadd_round_pd(V1(vC3), R, V1(vC2), RN);
    p = _mm512_fmadd_round_pd(R2, p, q23, RN);
    const __m512d q01 = _mm512_fmadd_round_pd(V1(vC1), R, V1(vC0), RN);
    const __m512d ShThR = _mm512_mul_round_pd(ShTh, R, RN);                            /* (Sh*Th)_h*R */
    p = _mm512_fmadd_round_pd(R2, p, q01, RN);
    lo = _mm512_fmadd_round_pd(p, ShThR, lo, RN);                                      /* Low + ShTh*R*Poly */
    const __m512d y = _mm512_add_round_pd(lo, ShTh, RN);
    const __m512d res = _mm512_scalef_round_pd(y, K, RN);                              /* * 2^K */
    return _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(res), sgn));
}

__attribute__((target("avx512f"))) void
svml_cbrt_ha_array(const double *x, double *y, size_t len)
{
    size_t i = 0;
    for (; i + 8 <= len; i += 8) {
        _mm512_storeu_pd(y + i, svml_cbrt8_ha(_mm512_loadu_pd(x + i)));
    }
    if (i < len) {
        const __mmask8 m = (__mmask8)((1u << (len - i)) - 1);
        _mm512_mask_storeu_pd(y + i, m, svml_cbrt8_ha(_mm512_maskz_loadu_pd(m, x + i)));
    }
}

#endif /* __x86_64__ */
