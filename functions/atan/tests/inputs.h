/* Test inputs for the double-precision atan tests: special values, every threshold of the vector
 * algorithm (table intervals j/4 +- 1/8, 7.875, 2^128) and of the scalar routine with neighbouring
 * doubles, and random values over [-8, 8], over log-uniform magnitudes and over all bit patterns. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng */

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 16384 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000123ULL), 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 1e-300, -1e-300, 5e-324, -5e-324,
        0x1p-1022, -0x1p-1022, 0x1.fffffffffffffp-1023, 1e308, -1e308, 0x1.fffffffffffffp+1023,
        -0x1.fffffffffffffp+1023, 0x1p128, 0x1p600, -0x1p600,
    };
    const double edges[] = {
        7.875, 3.875, 4.0, 4.125, 0x1p128,  /* vector: large-|x| path, second table half (j >= 16), Y clamp */
        0.125, 0x1p54, 0x1p-53, 0x1p-27, 0x1p-1022,  /* scalar routine */
        0x1p-26, 0x1p26, 1e10, 1e-10,
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
    for (int j = 0; j < 64; j++) {  /* table points j/4 and interval ends (2j+1)/8, +-3 ulps */
        for (int s = -3; s <= 3; s++) {
            buf[k++] = bitsd(dbits(j * 0.25) + (uint64_t)s * (j > 0));
            buf[k++] = -bitsd(dbits((2 * j + 1) * 0.125) + (uint64_t)s);
            buf[k++] = bitsd(dbits((2 * j + 1) * 0.125) + (uint64_t)s);
        }
    }
    for (size_t i = 0; i < n_random / 8 * 4; i++) {  /* uniform over [-8, 8] */
        buf[k++] = rng_uniform(-8.0, 8.0);
    }
    for (size_t i = 0; i < n_random / 8 * 2; i++) {  /* log-uniform magnitudes in [2^-60, 2^70], random sign */
        const uint64_t r = rng();
        double v = ldexp(1.0 + (double)(r >> 12) * 0x1p-52, (int)(r % 131) - 60);
        buf[k++] = (r & 2048) ? -v : v;
    }
    for (size_t i = 0; i < n_random / 8 * 2; i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
