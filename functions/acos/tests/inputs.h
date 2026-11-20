/* Test inputs for the double-precision acos tests: special values, every threshold of the vector
 * algorithm and of the scalar callout with neighbouring doubles, and random values over [-1, 1], just
 * outside it, over small magnitudes and over all bit patterns. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng */

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 8192 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000123ULL), 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 1e-300, -1e-300, 5e-324, -5e-324,
        0x1p-1022, -0x1p-1022, 0x1.fffffffffffffp-1023, 1e308, -1e308, 0x1.fffffffffffffp+1023,
        -0x1.fffffffffffffp+1023, 0.7071067811865476, -0.7071067811865476, 0.8660254037844386,
    };
    const double edges[] = {
        0.5, 1.0,                   /* vector: R = min(x^2, Y) switch, callout for |x| > 1 */
        0x1p-53, 0x1p-1022,         /* callout: tiny x */
        0x1p-26, 0x1p-27, 0.25, 0.75, 0.9375, 0.99, 0x1.fffffffp-1, 0x1p-1, 0x1.0000000000001p-1,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        for (int sg = 0; sg < 2; sg++) {
            const double v = sg ? -edges[i] : edges[i];
            for (int s = 0; s < 101; s++) {  /* the edge and 50 doubles on either side */
                buf[k++] = bitsd(dbits(v) + (uint64_t)(s - 50));
            }
        }
    }
    for (int j = 1; j <= 2000; j++) {  /* 1 - 2^-j*(1..)... : arguments approaching +-1 */
        const double v = 1.0 - ldexp(1.0 + (j % 7) * 0.125, -(j % 54) - 1);
        buf[k++] = (j & 1) ? v : -v;
    }
    for (size_t i = 0; i < n_random / 8 * 4; i++) {  /* uniform over [-1, 1] */
        buf[k++] = rng_uniform(-1.0, 1.0);
    }
    for (size_t i = 0; i < n_random / 8; i++) {  /* uniform over [-1.25, 1.25] */
        buf[k++] = rng_uniform(-1.25, 1.25);
    }
    for (size_t i = 0; i < n_random / 8; i++) {  /* log-uniform magnitudes in [2^-60, 1], random sign */
        const uint64_t r = rng();
        double v = ldexp(1.0 + (double)(r >> 12) * 0x1p-52, -(int)(r % 60) - 1);
        buf[k++] = (r & 2048) ? -v : v;
    }
    for (size_t i = 0; i < n_random / 8 * 2; i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
