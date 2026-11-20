/* Test inputs for the double-precision acosh tests: special values, the domain edge x = 1, every threshold
 * of the algorithm with neighbouring doubles, and random values over the interesting ranges and over all
 * bit patterns. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng, rng_uniform */

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 8192 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        1.0, -1.0, 0.5, 2.0, 1e-300, 5e-324, 0x1p-1022, 1e308, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
        0x1.6a09e667f3bccp+511, 0x1.6a09e667f3bcdp+511, 0x1p512, 0.75, 4.0 / 3.0, bitsd(0x3ff0000000000000ULL | 1),
        bitsd(0x3ff0000100000000ULL), bitsd(0x7ff0000100000000ULL), bitsd(0xbff0000000000000ULL),
    };
    const double edges[] = {
        1.0,                      /* domain edge: callout below, Yh = 0 at 1 */
        0x1.0003ffe002p+0,        /* ~sqrt(1 + 2^-11): small-input polynomial below */
        0x1p511,                  /* from here x/8 replaces x + sqrt(x^2 - 1) */
        0x1.fffffffffffffp+1023,  /* above: scalar callout (inf) */
        0x1.6a09e667f3bcdp+511,   /* x^2 overflows */
        0x1p26, 0x1p27,           /* x^2 - 1 rounds to x^2 */
        0x1.6a09e667f3bcdp+0, 2.0, 1.25,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) {
        buf[k++] = specials[i];
        buf[k++] = -specials[i];
    }
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        for (int s = 0; s < 201; s++) buf[k++] = bitsd(dbits(edges[i]) + (uint64_t)(s - 100));
    }
    for (int s = 0; s < 1000; s++) buf[k++] = bitsd(dbits(1.0) + (uint64_t)s * s * s);  /* just above 1 */
    /* Xin = x + sqrt(x^2 - 1) at the boundaries of the 4-bit reciprocal rounding and of the table:
     * Xin = 2^e * (1 + j/32), x = (Xin + 1/Xin)/2 */
    for (int e = 0; e <= 40; e += 3) {
        for (int j = 0; j < 32; j++) {
            const double xin = ldexp(1.0 + j / 32.0, e);
            const double x = 0.5 * (xin + 1.0 / xin);
            for (int s = -3; s <= 3; s++) buf[k++] = bitsd(dbits(x) + (uint64_t)s);
        }
    }
    const size_t q = n_random / 8;
    for (size_t i = 0; i < q; i++) buf[k++] = rng_uniform(1.0, 1.002);  /* around the small-input threshold */
    for (size_t i = 0; i < q; i++) buf[k++] = rng_uniform(0.5, 4.0);
    for (size_t i = 0; i < q; i++) buf[k++] = rng_uniform(1.0, 1000.0);
    for (size_t i = 0; i < q; i++) {  /* log-uniform positive magnitudes over the whole range */
        const uint64_t r = rng();
        buf[k++] = bitsd((r & 0x000fffffffffffffULL) | ((uint64_t)(r >> 52 & 0x7ff) % 0x7ff << 52));
    }
    for (size_t i = 0; i < n_random - 4 * q; i++) buf[k++] = bitsd(rng());  /* arbitrary bit patterns */
    return k;
}
