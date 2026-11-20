/* Test inputs for the double-precision tanh tests: special values, every threshold of the vector
 * algorithm (subinterval boundaries, saturation, callout) with neighbouring doubles, and random values
 * (uniform over the useful range, log-uniform magnitudes for the small-|x| region, arbitrary bit
 * patterns). Deterministic. */
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
        1.0, -1.0, 0.5, 1e-300, -1e-300, 5e-324, -5e-324, 0x1p-1022, -0x1p-1022, 0x1.fffffffffffffp-1023,
        0x1p-52, -0x1p-52, 0x1p-53, 0x1.fffffffffffffp-54, 0x1p-26, 1e308, -1e308,
        0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
    };
    const double edges[] = {
        0x1.8p+1023,              /* vector: callout for high word & 0x7ff80000 > 0x7fe00000 */
        0x1p+1023,
        19.061547465398498,       /* tanh rounds to 1 above ~atanh(1 - 2^-54) */
        18.714973875118524,       /* ~atanh(1 - 2^-53) */
        0x1p-26, 0x1p-54, 1e-300,
        0.5, 1.0, 2.0, 4.0, 8.0,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) {
        buf[k++] = specials[i];
    }
    for (int e = -3; e <= 4; e++) {  /* subinterval boundaries 1.5*2^e and 2^(e+1) (0.1875 ... 24) */
        const double b[2] = {ldexp(1.5, e), ldexp(1.0, e + 1)};
        for (int i = 0; i < 2; i++) {
            for (int s = 0; s < 101; s++) {
                const double v = bitsd(dbits(b[i]) + (uint64_t)(s - 50));
                buf[k++] = v;
                buf[k++] = -v;
            }
        }
    }
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        for (int s = 0; s < 101; s++) {  /* the edge and 50 doubles on either side, both signs */
            const double v = bitsd(dbits(edges[i]) + (uint64_t)(s - 50));
            buf[k++] = v;
            buf[k++] = -v;
        }
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [-30, 30] */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 60.0 - 30.0;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [-2, 2] */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 4.0 - 2.0;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* +-m*2^e, e uniform in [-60, 5], m uniform in [1, 2) */
        const uint64_t b = rng();
        const uint64_t e = (uint64_t)(1023 - 60) + (b >> 54) % 66;
        buf[k++] = bitsd((b << 63) | (e << 52) | (rng() & 0x000fffffffffffffULL));
    }
    for (size_t i = 0; i < n_random - 3 * (n_random / 4); i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
