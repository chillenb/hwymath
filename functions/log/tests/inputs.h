/* Test inputs for the double-precision log tests: special values, every threshold of the algorithm and of
 * the scalar callout with neighbouring doubles, and random values over the useful range and over all bit
 * patterns. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h"  /* dbits, bitsd, rng */

/* the 16 mantissas where the rounded vrcp14pd(m) changes (functions/log/hwy/rcp14_thresholds.h) */
static const uint64_t kRcpThresholdBits[16] = {
    0x3fff080000000000ULL, 0x3ffd41c000000000ULL, 0x3ffbad0000000000ULL, 0x3ffa41a000000000ULL,
    0x3ff8f9c000000000ULL, 0x3ff7d07000000000ULL, 0x3ff6c15000000000ULL, 0x3ff5c99000000000ULL,
    0x3ff4e5f000000000ULL, 0x3ff4143000000000ULL, 0x3ff3523000000000ULL, 0x3ff29e6000000000ULL,
    0x3ff1f70000000000ULL, 0x3ff15b4000000000ULL, 0x3ff0c97000000000ULL, 0x3ff040f000000000ULL,
};

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 40000 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000123ULL), bitsd(0xfff8000000000456ULL),
        1.0, -1.0, 0.5, 2.0, 10.0, 0.1, 1e-300, -1e-300, 5e-324, -5e-324, 0x1p-1022, 0x1.fffffffffffffp-1023,
        0x1p-1074, 1e308, -1e308, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023, 2.718281828459045,
        0x1p-1023, 0x1p+1023, 0x1p-537, 0x1.8p-1050,
    };
    const double edges[] = {
        1.0, 0.990234375, 1.009765625,  /* callout: |x - 1| = 0.009765625 boundary of the near-1 path */
        0.75, 1.5, 0x1p-1022, 0x1p+1023, 0x1.fffffffffffffp+1023, 0x1p-1074, 0x1p-1060, 0.5, 2.0,
        /* log2 callout: high-word windows around 1 */
        0.984375, 1.01953125, 0x1.ffff000000000p-1, 0x1.0000800000000p+0, 1.0 - 0x1p-31, 1.0 + 0x1p-31,
        1.0 - 0x1p-32, 1.0 + 0x1p-32, 1.0 - 0x1p-20, 1.0 + 0x1p-20,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        uint64_t b = dbits(edges[i]);
        for (int s = -50; s <= 50; s++) {  /* the edge and 50 doubles on either side */
            uint64_t v = b + (uint64_t)(int64_t)s;
            if ((v >> 63) == (b >> 63)) buf[k++] = bitsd(v);
        }
    }
    /* rounded-reciprocal thresholds and table boundaries (m = 1 + i/128) at several exponents */
    const int exps[] = {-1074 + 60, -1030, -1022, -60, -2, -1, 0, 1, 5, 100, 1023};
    for (size_t e = 0; e < sizeof exps / sizeof *exps; e++) {
        for (int j = 0; j < 16 + 129; j++) {
            const double m = j < 16 ? bitsd(kRcpThresholdBits[j]) : 1.0 + (j - 16) / 128.0;
            const uint64_t b = dbits(ldexp(m, exps[e]));
            for (int s = -8; s <= 8; s++) buf[k++] = bitsd(b + (uint64_t)(int64_t)s);
        }
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* uniform over [0, 4) */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 4.0;
    }
    for (size_t i = 0; i < n_random / 8; i++) {  /* uniform near 1: [0.98, 1.02) */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 0.04 + 0.98;
    }
    for (size_t i = 0; i < n_random / 4 - n_random / 8; i++) {  /* 1 +- 2^-[1, 60] * [1, 2) */
        const double m = 1.0 + (double)(rng() >> 12) * 0x1p-52;
        const uint64_t r = rng();
        buf[k++] = 1.0 + ldexp((r >> 63) ? -m : m, -1 - (int)(r % 60));
    }
    for (size_t i = 0; i < n_random / 4; i++) {  /* positive bit patterns: all exponents, incl. denormals */
        buf[k++] = bitsd(rng() >> 1);
    }
    for (size_t i = 0; i < n_random - 3 * (n_random / 4); i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
