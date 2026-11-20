/* Inputs for the vrcp14pd tests (test_rcp14.c: hardware/SDE vs the reference; hwy_test.cc: the Highway
 * emulation vs the reference): every 16-bit leading fraction pattern (vrcp14pd reads only 16 fraction
 * bits) at the exponents with special handling and some others, both signs, with and without random
 * low bits; denormals with every leading pattern; every exponent with random fractions; special values;
 * random bit patterns. buf must hold rcp14_max_inputs(n_random) values; returns how many. */
#include <stdint.h>

#include "hwymath/test_util.h"

static inline size_t
rcp14_max_inputs(size_t n_random)
{
    return 65536 * 22 + 2048 * 16 + 64 + n_random;
}

static size_t
rcp14_inputs(uint64_t *in, size_t n_random)
{
    static const uint64_t exps[] = {0x3ff, 0x3fe, 0x400, 0x001, 0x002, 0x7fc, 0x7fd, 0x7fe, 0x123, 0x5a5};
    size_t n = 0;
    for (size_t k = 0; k < sizeof exps / sizeof *exps; k++) {
        for (uint64_t f = 0; f < 65536; f++) {
            const uint64_t low = (k & 1) ? rng() & 0xfffffffffULL : 0; /* the 36 ignored bits */
            in[n++] = exps[k] << 52 | f << 36 | low;
            in[n++] = 0x8000000000000000ULL | exps[k] << 52 | f << 36 | (rng() & 0xfffffffffULL);
        }
    }
    for (uint64_t f = 0; f < 65536; f++) { /* denormals with every leading 16-bit pattern after normalization */
        const int sh = (int)(f % 52);
        in[n++] = (0x0010000000000000ULL | f << 36 | (rng() & 0xfffffffffULL)) >> (sh + 1);
        in[n++] = 0x8000000000000000ULL | (rng() >> (12 + sh));
    }
    for (uint64_t e = 0; e < 2048; e++) { /* every exponent */
        for (int j = 0; j < 8; j++) {
            in[n++] = e << 52 | (rng() >> 12);
            in[n++] = 0x8000000000000000ULL | e << 52 | (j ? rng() >> 12 : 0);
        }
    }
    static const uint64_t specials[] = {0, 0x8000000000000000ULL, 0x7ff0000000000000ULL, 0xfff0000000000000ULL,
                                        0x7ff8000000000000ULL, 0x7ff0000000000001ULL, 0xfff4000000000123ULL,
                                        0x0004000000000000ULL, 0x0004000000000001ULL, 0x8004000000000000ULL,
                                        0x0000000000000001ULL, 0x000fffffffffffffULL, 0x7fefffffffffffffULL,
                                        0x3ff0000000000000ULL, 0x7fe0000000000000ULL, 0x7fd0000000000000ULL,
                                        0x0008000000000000ULL, 0x0008000000000001ULL, 0x7fd0000000000001ULL};
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) in[n++] = specials[i];
    for (size_t i = 0; i < n_random; i++) in[n++] = rng();
    return n;
}
