/* Test inputs for the double-precision sinh tests: special values, every threshold of the vector
 * algorithm and of the scalar callout with neighbouring doubles, table boundaries, and random values
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
        0x1p-53,                  /* scalar: |x| < 2^-53 returns x */
        0x1.bdb8cdadbe120p-6,     /* scalar: odd Taylor polynomial below 0.0272 */
        0x1.57cd0e7026820p+4,     /* scalar: e^|x|/2 above 21.49 */
        0x1.633ce8fb9f87ep+9,     /* scalar: overflow above 710.476 */
        710.4758600739439,        /* log(2*DBL_MAX) */
        0x1.4ee3400000000p+9,     /* vector: callout for high word > 0x4084ee33 (669.78) */
        0.021660849392498291,     /* vector: N = 1/16 is replaced by N = 0 between ln2/32 */
        0.064982548177494873,     /*   and 3 ln2/32 */
        0.010830424696249145,     /* ln2/64 */
        0.69314718055994531, 1.0, 2.0, 20.0, 22.0, 650.0, 665.0, 669.0,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) {
        buf[k++] = specials[i];
    }
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        for (int s = 0; s < 101; s++) {  /* the edge and 50 doubles on either side, both signs */
            const double v = bitsd(dbits(edges[i]) + (uint64_t)(s - 50));
            buf[k++] = v;
            buf[k++] = -v;
        }
    }
    for (int j = 0; j < 64; j++) {  /* table boundaries of both algorithms: (M*64 + j + 1/2)*ln2/64 */
        for (int m = 0; m <= 1024; m += 37) {
            const double v = (m * 64 + j + 0.5) * 0.010830424696249145;
            buf[k++] = v;
            buf[k++] = -bitsd(dbits(v) + 1);
        }
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [-720, 720] */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 1440.0 - 720.0;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [-25, 25] */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 50.0 - 25.0;
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* +-m*2^e, e uniform in [-60, 3], m uniform in [1, 2) */
        const uint64_t b = rng();
        const uint64_t e = (uint64_t)(1023 - 60) + (b >> 54) % 64;
        buf[k++] = bitsd((b << 63) | (e << 52) | (rng() & 0x000fffffffffffffULL));
    }
    for (size_t i = 0; i < n_random - 3 * (n_random / 4); i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
