/* Test inputs for the double-precision atanh tests: special values, the domain edges |x| = 1, every
 * threshold of the algorithm with neighbouring doubles, and random values over the interesting ranges and
 * over all bit patterns. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng, rng_uniform */

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 16384 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, INFINITY, NAN, bitsd(0x7ff0000000000001ULL), bitsd(0x7ff4000000000000ULL), 1.0, 2.0, 0.5, 1e-300,
        5e-324, 0x1p-1022, 0x1.fffffffffffffp-1023, 0x1p-26, 0x1p-27, 0x1p-53, 0x1p-54, 0x1.fffffffffffffp-1,
        0x1.ffffffffffffep-1, 1e308, 0x1.fffffffffffffp+1023, bitsd(0x3ff0000000000001ULL), bitsd(0x3ff0000100000000ULL),
        bitsd(0x7ff0000100000000ULL), 0.75, 0.25, 1.0 / 3.0,
    };
    const double edges[] = {
        1.0,   /* domain edge: callout from here */
        0.5,   /* 1 - |x| changes binade */
        0x1p-53, 0x1p-54,  /* 1 - |x| and 1 + |x| round to 1 */
        0.25, 0.75, 0x1p-5,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) {
        buf[k++] = specials[i];
        buf[k++] = -specials[i];
    }
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        for (int s = 0; s < 201; s++) {
            const double v = bitsd(dbits(edges[i]) + (uint64_t)(s - 100));
            buf[k++] = v;
            buf[k++] = -v;
        }
    }
    /* boundaries of the 4-bit reciprocal rounding and of the tables: 1 + |x| = 1 + j/32, 1 - |x| = 2^-e (1 + j/32) */
    for (int j = 0; j < 32; j++) {
        for (int s = -3; s <= 3; s++) buf[k++] = bitsd(dbits(j / 32.0) + (uint64_t)s);
        for (int e = 1; e <= 52; e += 3) {
            const double x = 1.0 - ldexp(1.0 + j / 32.0, -e);
            for (int s = -3; s <= 3; s++) buf[k++] = bitsd(dbits(x) + (uint64_t)s);
        }
    }
    const size_t q = n_random / 8;
    for (size_t i = 0; i < 2 * q; i++) buf[k++] = rng_uniform(-1.0, 1.0);
    for (size_t i = 0; i < q; i++) {  /* close to +-1: 1 - |x| log-uniform in [2^-53, 1] (built from bits: no libm) */
        const uint64_t r = rng();
        const double v = 1.0 - bitsd(((uint64_t)(1023 - 1 - (r >> 58) % 53) << 52) | (r & 0x000fffffffffffffULL));
        buf[k++] = (r >> 57 & 1) ? v : -v;
    }
    for (size_t i = 0; i < q; i++) {  /* log-uniform magnitudes below 1, random sign */
        const uint64_t r = rng();
        buf[k++] = bitsd((r & 0x800fffffffffffffULL) | ((uint64_t)(r >> 52 & 0x7ff) % 0x3ff << 52));
    }
    for (size_t i = 0; i < n_random - 4 * q; i++) buf[k++] = bitsd(rng());  /* arbitrary bit patterns */
    return k;
}
