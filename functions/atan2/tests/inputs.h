/* Test inputs for the double-precision atan2 tests, as pairs (y[i], x[i]): all combinations of special
 * values, ratios |y/x| around every interval threshold of the vector algorithm, magnitudes around the
 * vector path's exponent limits (2^-1020, 2^993), and random pairs (uniform, log-uniform magnitudes,
 * random ratios, arbitrary bit patterns). */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng */

#define INPUTS_EXTRA 16384  /* make_inputs writes at most n_random + INPUTS_EXTRA pairs */

static size_t
make_inputs(double *ys, double *xs, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000123ULL), 1.0, -1.0, 0.5, -2.0, 5e-324, -5e-324, 0x1p-1030, -0x1.8p-1040, 0x1p-1022,
        -0x1p-1022, 0x1p-1020, 0x1.fffffffffffffp-1021, 0x1p993, -0x1.fffffffffffffp992, 0x1.fffffffffffffp+1023,
        -0x1.fffffffffffffp+1023, 1e300, -3.0, 0.4375,
    };
    const size_t ns = sizeof specials / sizeof *specials;
    for (size_t i = 0; i < ns; i++) {
        for (size_t j = 0; j < ns; j++) { ys[k] = specials[i]; xs[k] = specials[j]; k++; }
    }
    /* |y/x| around the interval thresholds 7/16, 11/16, 19/16, 39/16 (the comparisons are ax*t < ay) */
    const double thr[] = {0.4375, 0.6875, 1.1875, 2.4375};
    const double xv[] = {1.0, 3.3, 0x1.fffffffffffffp-1, 1.2345e-5, 7.77e5, 0x1p-1000, 0x1p990};
    for (size_t t = 0; t < 4; t++) {
        for (size_t v = 0; v < sizeof xv / sizeof *xv; v++) {
            const double yc = xv[v] * thr[t];
            for (int s = -20; s <= 20; s++) {
                const int sg = (s + 40) & 3;
                ys[k] = (sg & 1 ? -1 : 1) * bitsd(dbits(yc) + (uint64_t)(int64_t)s);
                xs[k] = (sg & 2 ? -1 : 1) * xv[v];
                k++;
            }
        }
    }
    /* magnitudes around the exponent limits of the vector path, in ulps and in steps of the high word */
    const double lim[] = {0x1p-1020, 0x1p993};
    const double other[] = {1.0, 0x1p-1020, 0x1p993, 0x1p-1000, 0x1p980, 3.0, -0.75};
    for (size_t l = 0; l < 2; l++) {
        for (int s = -30; s <= 30; s++) {
            const double v = bitsd(dbits(lim[l]) + (uint64_t)(int64_t)s * ((s & 1) ? 1 : 0x100000000ULL));
            for (size_t o = 0; o < sizeof other / sizeof *other; o++) {
                ys[k] = (s & 2) ? -v : v; xs[k] = other[o]; k++;
                ys[k] = other[o]; xs[k] = (s & 4) ? -v : v; k++;
            }
        }
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [-2, 2]^2 */
        ys[k] = rng_uniform(-2.0, 2.0); xs[k] = rng_uniform(-2.0, 2.0); k++;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* y = x * ratio, ratio log-uniform in [2^-8, 2^8] */
        /* (one rng() call per statement: the evaluation order of function arguments is unspecified) */
        const double x = rng_uniform(-100.0, 100.0);
        const double m = 1.0 + (double)(rng() >> 12) * 0x1p-52;
        const double q = ldexp(m, (int)(rng() % 16) - 8);
        ys[k] = (rng() & 1) ? -x * q : x * q; xs[k] = x; k++;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* log-uniform magnitudes in [2^-1080, 2^1024), random signs */
        const uint64_t r = rng();
        const int ey = (int)(rng() % 2104) - 1080;
        const double mx = 1.0 + (double)(rng() >> 12) * 0x1p-52;
        const int ex = (int)(rng() % 2104) - 1080;
        ys[k] = ldexp(1.0 + (double)(r >> 12) * 0x1p-52, ey) * ((r & 1) ? -1 : 1);
        xs[k] = ldexp(mx, ex) * ((r & 2) ? -1 : 1);
        k++;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* arbitrary bit patterns */
        ys[k] = bitsd(rng()); xs[k] = bitsd(rng()); k++;
    }
    return k;
}
