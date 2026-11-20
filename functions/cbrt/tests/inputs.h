/* Test inputs for the double-precision cbrt tests: special values, powers of two over the whole exponent
 * range, perfect cubes, denormals, every threshold of the algorithm (the rounded-reciprocal thresholds and
 * the exact boundaries 32/(2k+1), table boundaries of the scalar routine) with neighbouring doubles, and
 * random values. Deterministic. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng */
#include "../rcp14_model.h"

#define CBRT_FIXED_INPUTS 32768  /* upper bound on the non-random inputs */

/* adds v and its neighbours v +- 1..w ulps (same sign) */
static size_t
add_around(double *buf, size_t k, double v, int w)
{
    for (int s = -w; s <= w; s++) buf[k++] = bitsd(dbits(v) + (uint64_t)(int64_t)s);
    return k;
}

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + CBRT_FIXED_INPUTS. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN,
        bitsd(0x7ff0000000000001ULL), bitsd(0xfff0000000000001ULL),  /* signalling NaNs */
        bitsd(0x7ff4000000000005ULL), bitsd(0x7ff8000000000003ULL),  /* NaN payloads with various low bits */
        bitsd(0xfff8000000000006ULL), bitsd(0x7fffffffffffffffULL), bitsd(0xffffffffffffffffULL),
        1.0, -1.0, 8.0, -8.0, 27.0, 0.125, 1e-300, -1e-300, 1e300, 5e-324, -5e-324, 1e-310, -1e-310,
        0x1p-1022, -0x1p-1022, 0x1.fffffffffffffp-1023, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
        0x1p-1074, 0x1p-1073, 0x1.8p-1073, 3.0, 2.0, 4.0, 0.1, -0.1, 1e-5, 1e5, 2.5, 1000.0, 1e-3,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];

    /* powers of two over the whole range (including denormals) and their neighbours, both signs */
    for (int e = -1074; e <= 1023; e++) {
        const double p = ldexp(1.0, e);
        k = add_around(buf, k, p, 1);
        buf[k++] = -p;
    }
    /* perfect cubes n^3 (exact), scaled by 2^(3s): normal, tiny normal, denormal, huge */
    const int scales[] = {0, -330, -355, 300};
    for (int n = 1; n <= 1536; n++) {
        const double c = (double)n * n * n;
        for (int s = 0; s < 4; s++) {
            const double v = ldexp(c, 3 * scales[s]);
            buf[k++] = (n & 1) ? -v : v;
        }
    }
    /* reciprocal thresholds: Intel's (multiples of 2^-16) and the exact 32/(2k+1), with exponents
     * covering E mod 3 = 0, 1, 2 and a denormal */
    const int exps[] = {0, 1, 2, -3, 31, -1040};
    for (int i = 0; i < 8; i++) {
        const double b[2] = {cbrt_rcp14_threshold[i], 32.0 / (2 * (15 - i) + 1)};
        for (int j = 0; j < 2; j++) {
            for (size_t e = 0; e < sizeof exps / sizeof *exps; e++) {
                k = add_around(buf, k, ldexp(b[j], exps[e]), 20);
            }
        }
        k = add_around(buf, k, cbrt_rcp14_threshold[i] - 0x1p-16, 2);  /* start of the previous prefix */
    }
    /* the scalar routine's table boundaries 1 + i/32 and the mantissa boundaries 1 and 2 */
    for (int i = 0; i <= 32; i++) {
        k = add_around(buf, k, 1.0 + i / 32.0, 5);
        k = add_around(buf, k, ldexp(1.0 + i / 32.0, -1050), 5);
    }
    /* random: half arbitrary bit patterns, a quarter with uniform mantissa and E in [-6, 6), a quarter
     * uniform in [-1000, 1000] */
    for (size_t i = 0; i < n_random / 2; i++) buf[k++] = bitsd(rng());
    for (size_t i = 0; i < n_random / 4; i++) {
        const uint64_t r = rng();
        buf[k++] = bitsd((r & 0x800fffffffffffffULL) | ((uint64_t)(1023 - 6 + (int)((r >> 52) % 12)) << 52));
    }
    for (size_t i = n_random / 2 + n_random / 4; i < n_random; i++) buf[k++] = rng_uniform(-1000.0, 1000.0);
    return k;
}
