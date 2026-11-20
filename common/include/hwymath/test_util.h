/* Helpers shared by the tests of all functions: bit casts between double and uint64_t, and a small
 * deterministic random number generator (identical sequences on every platform). C and C++. */
#ifndef HWYMATH_TEST_UTIL_H
#define HWYMATH_TEST_UTIL_H

#include <stdint.h>
#include <string.h>

static inline uint64_t
dbits(double d)
{
    uint64_t u;
    memcpy(&u, &d, sizeof u);
    return u;
}

static inline double
bitsd(uint64_t u)
{
    double d;
    memcpy(&d, &u, sizeof d);
    return d;
}

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;

/* xorshift64* */
static inline uint64_t
rng(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545f4914f6cdd1dULL;
}

/* uniform double in [lo, hi) from the generator above (compile with -ffp-contract=off for identical
 * results on every architecture) */
static inline double
rng_uniform(double lo, double hi)
{
    return ((double)(rng() >> 11) * 0x1p-53) * (hi - lo) + lo;
}

#endif
