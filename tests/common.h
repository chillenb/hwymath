/* Test inputs shared by the exp tests: special values, every threshold of the algorithm with
 * neighbouring doubles, and random values over the interesting range and over all bit patterns. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
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
static inline uint64_t
rng(void)  /* xorshift64* */
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545f4914f6cdd1dULL;
}

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 4096 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        1.0, -1.0, 0.5, 1e-300, -1e-300, 5e-324, -5e-324, 0x1p-52, -0x1p-52, 0x1p-53, 0x1.fffffffffffffp-53,
        0x1p-1022, 1e308, -1e308, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
    };
    const double edges[] = {
        709.782712893384,         /* scalar overflow threshold log(DBL_MAX) */
        707.7032713517042,        /* vector callout threshold */
        -707.7032713517042,
        -708.3964185322641,       /* scalar: denormal results below log(DBL_MIN) */
        -745.1332191019411,       /* scalar: result rounds to zero below this */
        -744.4400719213812,       /* log(smallest denormal) */
        709.0895657128241,        /* ~1023*ln2: M = 1023/1024 boundary */
        0.6931471805599453, -0.6931471805599453, 0.04332169878499658, 0.0, 1.0,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        double v = edges[i];
        for (int s = 0; s < 101; s++) {  /* the edge and 50 doubles on either side */
            buf[k++] = bitsd(dbits(v) + (uint64_t)(s - 50));
        }
    }
    for (int j = 0; j < 64; j++) {  /* table boundaries of both algorithms: (M*64 + j)*ln2/64 */
        for (int m = -1080; m <= 1024; m += 131) {
            buf[k++] = (m * 64 + j) * 0.010830424696249145;
        }
    }
    for (size_t i = 0; i < n_random / 2; i++) {  /* uniform over [-750, 750] */
        buf[k++] = ((double)(rng() >> 11) * 0x1p-53) * 1500.0 - 750.0;
    }
    for (size_t i = 0; i < n_random / 2; i++) {  /* arbitrary bit patterns */
        buf[k++] = bitsd(rng());
    }
    return k;
}
