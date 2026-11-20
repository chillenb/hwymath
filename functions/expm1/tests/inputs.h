/* Test inputs for the double-precision expm1 tests: special values, every threshold of the algorithm with
 * neighbouring doubles, table boundaries, and random values over the interesting range and over all bit
 * patterns. Deterministic (same inputs on every platform). */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng, rng_uniform */

#define INPUTS_EXTRA 32768  /* non-random inputs: buf must hold n_random + INPUTS_EXTRA values */

/* Fills buf with test inputs; returns how many. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000123ULL), bitsd(0xfff0000000000001ULL),
        1.0, -1.0, 0.5, -0.5, 2.0, 10.0, -10.0, 1e-300, -1e-300, 5e-324, -5e-324, 0x1p-1022, -0x1p-1022,
        0x1p-52, -0x1p-52, 0x1p-53, -0x1p-53, 1e308, -1e308, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
        -1e18, -1e100, -0x1p1023, 0x1p1023,
    };
    const double ln2 = 0.6931471805599453;
    const double edges[] = {
        707.7032713517042,        /* vector callout threshold (x >= this) */
        709.782712893384,         /* scalar: overflow threshold log(DBL_MAX) */
        709.0895657128241,        /* scalar: M = 1023/1024 boundary */
        -128 * ln2, -128.03125 * ln2, -127.96875 * ln2,  /* N = -128: clamp of the scaling */
        -37.42994775023705,       /* expm1(x) rounds to -1 below about this */
        -36.7368005696771,        /* ~ -53 ln2 */
        -ln2 / 32, ln2 / 32,      /* N = 0 / -1/16 boundary: Th' < 1 switches the Fast2Sum order */
        -ln2 / 16, ln2 / 16, -ln2, ln2, -1.0, 1.0,
        0.0,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        double v = edges[i];
        for (int s = 0; s < 101; s++) {  /* the edge and 50 doubles on either side (crossing 0 for 0) */
            const int64_t o = s - 50;
            buf[k++] = v == 0.0 ? (o < 0 ? -bitsd((uint64_t)-o) : bitsd((uint64_t)o)) : bitsd(dbits(v) + (uint64_t)o);
        }
    }
    for (int e = -1074; e <= 10; e++) {  /* powers of two of both signs */
        buf[k++] = -ldexp(1.0, e);
        buf[k++] = ldexp(1.0, e);
        buf[k++] = -ldexp(1.0 + 0x1p-52, e);
        buf[k++] = ldexp(1.0 - 0x1p-53, e);
    }
    for (int j = 0; j < 32; j++) {  /* table and rounding boundaries (M*16 + j/2)*ln2/16, and neighbours */
        for (int m = -1080; m <= 1024; m += 23) {
            const double v = (m * 16 + j * 0.5) * (ln2 / 16);
            buf[k++] = v;
            buf[k++] = nextafter(v, -INFINITY);
            buf[k++] = nextafter(v, INFINITY);
        }
    }
    for (size_t i = 0; i < n_random / 4; i++) buf[k++] = rng_uniform(-750.0, 750.0);
    for (size_t i = 0; i < n_random / 4; i++) buf[k++] = rng_uniform(-1.0, 1.0);
    for (size_t i = 0; i < n_random / 2; i++) buf[k++] = bitsd(rng());  /* arbitrary bit patterns */
    return k;
}
