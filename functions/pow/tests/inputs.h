/* Test inputs (x, y) for the double-precision pow tests: grids of special values, every threshold of the
 * vector algorithm and of the scalar callout with neighbouring doubles, the reciprocal switch points, |x|
 * near 1 with huge |y|, random pairs over the useful ranges and random bit patterns. Deterministic. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng, rng_uniform */
#include "../hwy/rcp14_switch.inc"

#define POW_INPUTS_FIXED 200000 /* upper bound of the non-random inputs */

/* nearest double to v moved by k ulps (k may be negative); stays within finite doubles of v's sign */
static inline double
ulps_away(double v, int k)
{
    if (v == 0.0) return k >= 0 ? bitsd((uint64_t)k) : -bitsd((uint64_t)-k);
    const uint64_t b = dbits(v);
    return bitsd(b + (uint64_t)(int64_t)(v > 0 ? k : -k));
}

/* log2 and 2^u from + - * / only, so that the inputs are bit-identical on every architecture (libm's log2
 * and exp2 may differ); accurate to ~1e-15, enough to place y within a few ulps of a threshold */
static inline double
det_log2(double v)
{
    int e;
    double m = frexp(v, &e);  /* [1/2, 1) */
    if (m < 0.70710678118654752) { m *= 2; e--; }
    const double s = (m - 1) / (m + 1), s2 = s * s;
    double t = 0, p = s;
    for (int i = 1; i < 40; i += 2) { t += p / i; p *= s2; }
    return e + 2 * t * 1.4426950408889634;
}
static inline double
det_exp2(double u)  /* roughly 2^u (log-uniform when u is uniform) */
{
    const double e = floor(u);
    return ldexp(1.0 + (u - e), (int)e);
}

/* Fills x[], y[] with test pairs; returns how many. The arrays must hold n_random + POW_INPUTS_FIXED. */
static size_t
make_inputs(double *x, double *y, size_t n_random)
{
    size_t k = 0;
#define ADD(a, b) do { x[k] = (a); y[k] = (b); k++; } while (0)
    const double sx[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 3.0, -3.0, 10.0, -10.0, 0.1, M_E, 1.5, 0.75,
        0x1.0000000000001p0, 0x1.fffffffffffffp-1, -0x1.0000000000001p0, -0x1.fffffffffffffp-1,
        5e-324, -5e-324, 0x1p-1022, -0x1p-1022, 0x0.fffffffffffffp-1022, 0x1p-1060, 1e-300, -1e-300,
        1e300, -1e300, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023, 0x1p-1074 * 3, 7.0, -7.0,
    };
    const double sy[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5, 1.5, -1.5, 0.25, 4.0, -4.0, 1023.0, 1024.0, -1074.0, -1075.0,
        0x1p20 + 1, 0x1p21 + 1, 0x1p21 + 2, 0x1p31 + 1, 0x1p32 + 1, 0x1p52 - 1, 0x1p52 + 1, 0x1p53, 0x1p53 + 2,
        -(0x1p52 + 1), 0x1p52 - 0.5, 0x1p51 + 0.5, 0x1p20 + 0.5, 0x1p60, 1e300, -1e300,
        5e-324, -5e-324, 0x1p-1022, 1e-300, 1e-20, -1e-20, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
    };
    const size_t nsx = sizeof sx / sizeof *sx, nsy = sizeof sy / sizeof *sy;
    for (size_t i = 0; i < nsx; i++) {
        for (size_t j = 0; j < nsy; j++) ADD(sx[i], sy[j]);
    }
    /* integer y of every exponent (odd and even, positive and negative) with negative and positive x */
    for (int e = 0; e <= 60; e++) {
        const double b = ldexp(1.0, e);
        const double ys[] = {b, b + 1, b - 1, b + 2, b + 0.5, -b, -(b + 1), -(b + 0.5)};
        for (size_t j = 0; j < sizeof ys / sizeof *ys; j++) {
            ADD(-1.0000001, ys[j]); ADD(-0.9999999, ys[j]); ADD(-1.0, ys[j]); ADD(-2.0, ys[j]); ADD(1.0000001, ys[j]);
        }
    }
    /* thresholds in z = y*log2(x): vector callout (|S| > 1021.5), overflow (1024), denormal results
     * (-1022 .. -1032 .. -1084 in the scalar routine), underflow (-1074, -1075), exp table boundaries, and the
     * scalar routine's tiny/huge shortcuts (exponent of log2|x| + exponent of y <= 1984, >= 2057) */
    const double xs[] = {2.0, 0.5, M_E, 10.0, 1.5, 0.75, 1.0 + 0x1p-20, 1.0 - 0x1p-30, 1e-300, 1e300, 0x1p-1060,
                         0x1.fffffffffffffp+1023, 1.0 + 0x1p-52, 7.0, -2.0, -0.5, -10.0};
    const double zs[] = {1021.5, -1021.5, 1024.0, -1022.0, -1032.0, -1074.0, -1075.0, -1084.0, -1085.0, 1023.0,
                         1023.5, -1022.5, 0.0625, -0.0625, 0.5, 1.0, 0x1p-62, -0x1p-62, 0x1p-61, 0x1p-63, 0x1p11,
                         -0x1p11, 0x1p12};
    for (size_t i = 0; i < sizeof xs / sizeof *xs; i++) {
        const double l2 = det_log2(fabs(xs[i]));
        for (size_t j = 0; j < sizeof zs / sizeof *zs; j++) {
            double yv = zs[j] / l2;
            if (xs[i] < 0) yv = nearbyint(yv) + (j & 1);  /* integer y for negative x */
            for (int s = -40; s <= 40; s++) ADD(xs[i], ulps_away(yv, s));
            if (xs[i] < 0) for (int s = -3; s <= 3; s++) ADD(xs[i], yv + s);
        }
    }
    /* powers of two and exp table boundaries: z = y exactly for x = 2 */
    for (int j = -16 * 1100; j <= 16 * 1100; j += 7) {
        ADD(2.0, j / 16.0);
        ADD(2.0, ulps_away(j / 16.0, 1));
        ADD(2.0, ulps_away(j / 16.0, -1));
    }
    /* x near 1 (few ulps) with large |y|, and x = 1 +- tiny with y moderate */
    for (int s = -20; s <= 20; s++) {
        if (s == 0) continue;
        const double xv = ulps_away(1.0, s);
        for (int e = 0; e <= 62; e += 2) {
            const double yv = ldexp(1.0, e) * (1.0 + 0.37 * (e & 3));
            ADD(xv, yv); ADD(xv, -yv); ADD(xv, ldexp(1.0, e));
        }
    }
    /* reciprocal switch points (rcp14 rounding) and table boundaries of the 32-entry log table:
     * mantissas at the switch prefixes +-3 with low bits all 0, all 1 and random, at several exponents */
    static const unsigned sw[32] = {POW_RCP14_SWITCH_PREFIXES};
    for (int i = 0; i < 32; i++) {
        for (int dp = -3; dp <= 3; dp++) {
            const uint64_t pre = (uint64_t)(sw[i] + dp) & 0xffff;
            const uint64_t lows[3] = {0, (1ULL << 36) - 1, rng() & ((1ULL << 36) - 1)};
            for (int l = 0; l < 3; l++) {
                const double m = bitsd(0x3fe0000000000000ULL | pre << 36 | lows[l]);
                const int es[] = {0, 1, 2, -1, 5, -30, 600, -1000};
                for (size_t e = 0; e < sizeof es / sizeof *es; e++) {
                    const double xv = ldexp(m, es[e]);
                    ADD(xv, 0.7); ADD(xv, -13.25); ADD(xv, 1.0); ADD(xv, rng_uniform(-30, 30));
                }
                ADD(ldexp(m, -1050), 0.3);  /* denormal x */
            }
        }
    }
    for (int i = 0; i <= 64; i++) {  /* m = k/64 boundaries: rcp exactly 1/m */
        for (int s = -2; s <= 2; s++) {
            const double m = ulps_away((64.0 + i) / 128.0, s);
            ADD(m, 2.5); ADD(2 * m, -0.3); ADD(ldexp(m, 100), 1.7);
        }
    }
    /* scalar callout: its 3-level reciprocal table boundaries (x mantissa bits 14-19, 10-15, 4-11) */
    for (int i = 0; i < 4096; i++) {
        const double xv = bitsd(0x3ff0000000000000ULL | (uint64_t)i << 40 | (rng() & 0xffffffffffULL));
        ADD(-xv, 3.0); ADD(-xv, -5.0); ADD(xv, 800.0); ADD(xv * 0x1p-1060, -0.9);
    }
#undef ADD
    /* random pairs */
    const size_t nr = n_random / 8;
    for (size_t i = 0; i < nr; i++) {  /* x = 2^u over all positive doubles, z uniform in [-1100, 1100] */
        const double xv = det_exp2(rng_uniform(-1074, 1024));
        if (xv == 0.0 || isinf(xv) || xv == 1.0) continue;
        x[k] = xv; y[k] = rng_uniform(-1100, 1100) / det_log2(xv); k++;
    }
    for (size_t i = 0; i < nr; i++) {  /* x in (0, 4), y in [-100, 100] */
        x[k] = rng_uniform(0, 4); y[k] = rng_uniform(-100, 100); k++;
    }
    for (size_t i = 0; i < nr; i++) {  /* |x| near 1, y large */
        x[k] = 1.0 + rng_uniform(-0x1p-20, 0x1p-20);
        const double m = rng_uniform(-1, 1);  /* separate statements: argument evaluation order is unspecified */
        y[k] = ldexp(m, (int)(rng() % 60)); k++;
    }
    for (size_t i = 0; i < nr; i++) {  /* z small: |y*log2 x| < 2^-40 down to tiny */
        x[k] = rng_uniform(0.01, 100);
        const double m = rng_uniform(-1, 1);
        y[k] = ldexp(m, -(int)(rng() % 1100)); k++;
    }
    for (size_t i = 0; i < nr; i++) {  /* negative x with integer y, and some non-integer y */
        const double xv = -det_exp2(rng_uniform(-60, 60));
        double yv = nearbyint(rng_uniform(-1000, 1000) / fabs(det_log2(-xv)));
        if (i % 8 == 0) yv += 0.5;
        x[k] = xv; y[k] = yv; k++;
    }
    for (size_t i = 0; i < nr; i++) {  /* z near the overflow / underflow / denormal boundaries */
        const double xv = det_exp2(rng_uniform(-1074, 1024));
        if (xv == 1.0 || xv == 0.0 || isinf(xv)) continue;
        x[k] = xv; y[k] = rng_uniform(-1090, -1015) * (i & 1 ? -1 : 1) / det_log2(xv); k++;
    }
    for (size_t i = 0; i < n_random - 6 * nr; i++) {  /* arbitrary bit patterns */
        x[k] = bitsd(rng());
        y[k] = (i & 1) ? bitsd(rng()) : bitsd(rng() & 0xc07fffffffffffffULL);  /* half of y in [2^-1023, 2^9) */
        k++;
    }
    return k;
}
