/* Test inputs for the double-precision sin, cos and tan tests (shared by functions/{sin,cos,tan}/tests):
 * special values, every threshold of the algorithms with neighbouring doubles, the rounding boundaries
 * of the main-path reduction, doubles close to multiples of pi/2 and pi/16, known hard cases for the
 * large-argument reduction, and random values (uniform over several ranges, log-uniform over all
 * exponents, arbitrary bit patterns). Deterministic. */
#include <math.h>
#include <stdlib.h>

#include "hwymath/test_util.h" /* dbits, bitsd, rng */

/* pi/2 = PIO2_1 + PIO2_2 + PIO2_3 (fdlibm's 33+33+53-bit split) */
#define TRIG_PIO2_1 0x1.921fb544p+0
#define TRIG_PIO2_2 0x1.0b4611a6p-34
#define TRIG_PIO2_3 0x1.3198a2e037073p-69

/* the double next to k*pi/2 (k*PIO2_1 and k*PIO2_2 are exact for k < 2^20; close enough above) */
static inline double
trig_near_k_pio2(double k)
{
    return k * TRIG_PIO2_1 + (k * TRIG_PIO2_2 + k * TRIG_PIO2_3);
}

/* Fills buf with test inputs; returns how many. buf must hold at least n_random + 40000 values. */
static size_t
make_inputs(double *buf, size_t n_random)
{
    size_t k = 0;
    const double specials[] = {
        0.0, -0.0, INFINITY, -INFINITY, NAN, -NAN, bitsd(0x7ff0000000000001ULL), bitsd(0xfff4000000000000ULL),
        bitsd(0x7ff8000000000001ULL), 1.0, -1.0, 0.5, 2.0, 3.0, 10.0, 100.0, 1e-300, -1e-300, 5e-324, -5e-324,
        0x1p-1022, -0x1p-1022, 0x1.fffffffffffffp-1023, 0x1p-52, 0x1p-26, 0x1p-20, 0x1p-10, 1e22, 1e100, 1e200,
        1e300, 1e308, -1e308, 0x1.fffffffffffffp+1023, -0x1.fffffffffffffp+1023,
        0x1.6ac5b262ca1ffp+849,  /* 6381956970095103*2^797: the double closest to a multiple of pi/2 */
        -0x1.6ac5b262ca1ffp+849,
        0x1.921fb54442d18p+0, 0x1.921fb54442d18p+1, 0x1.921fb54442d18p+2,  /* pi/2, pi, 2pi */
        0x1.921fb54442d18p-1, 0x1.921fb54442d18p-4,                         /* pi/4, pi/32 */
    };
    const double edges[] = {
        0x1p+24, -0x1p+24,                      /* sin/cos: large-argument threshold */
        0x1.000000e4db24cp+16, -0x1.000000e4db24cp+16, /* tan: large-argument threshold */
        0x1p-20,                                /* large path: smaller |x| pass through */
        0x1p+16, 0x1p+52, 0x1p+53, 0x1p+63, 0x1p+64, 0x1p+1023,
        0x1.6ac5b262ca1ffp+849, 0x1p-1022, 0x1p-1074 * 50,
    };
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) buf[k++] = specials[i];
    for (size_t i = 0; i < sizeof edges / sizeof *edges; i++) {
        const double v = edges[i];
        for (int s = 0; s < 101; s++) buf[k++] = bitsd(dbits(v) + (uint64_t)(s - 50));  /* +-50 doubles */
    }
    /* the doubles next to k*pi/2 (zeros and poles), +-2 ulps, for k up to 2^20 and some larger k */
    for (int m = -2000; m <= 2000; m++) {
        const double v = trig_near_k_pio2((double)m);
        for (int s = -2; s <= 2; s++) buf[k++] = bitsd(dbits(v) + (uint64_t)(int64_t)s);
    }
    for (int i = 0; i < 2000; i++) {
        const double kk = (double)(rng() >> (64 - 20 - (i % 30)));  /* up to 2^49 */
        const double v = trig_near_k_pio2(kk);
        for (int s = -1; s <= 1; s++) buf[k++] = bitsd(dbits(v) + (uint64_t)(int64_t)s);
        buf[k++] = -v;
    }
    /* rounding boundaries of the main-path reduction: N = round(x*16/pi) changes at (n + 1/2)*pi/16 */
    for (int n = -200; n <= 200; n++) {
        const double v = trig_near_k_pio2((n + 0.5) / 8.0);
        for (int s = -2; s <= 2; s++) buf[k++] = bitsd(dbits(v) + (uint64_t)(int64_t)s);
    }
    /* rounding boundaries of the large-path reduction, (2M + 1)*pi/512, for |x| ~ 2^17 (tan) and ~ 2^25 */
    for (int n = 0; n < 1024; n++) {
        const double v1 = trig_near_k_pio2((0x1p25 + 2 * n + 1) / 256.0);
        const double v2 = trig_near_k_pio2((0x1p33 + 2 * n + 1) / 256.0);
        buf[k++] = v1;
        buf[k++] = bitsd(dbits(v1) + 1);
        buf[k++] = v2;
        buf[k++] = bitsd(dbits(v2) + 1);
    }
    for (size_t i = 0; i < n_random / 8; i++) buf[k++] = rng_uniform(-10.0, 10.0);
    for (size_t i = 0; i < n_random / 8; i++) buf[k++] = rng_uniform(-1e5, 1e5);
    for (size_t i = 0; i < n_random / 8; i++) buf[k++] = rng_uniform(-0x1p24, 0x1p24);
    for (size_t i = 0; i < n_random / 8; i++) buf[k++] = rng_uniform(-0x1p30, 0x1p30);
    /* log-uniform: random mantissa and sign, exponent uniform in [-1074, 1023] (biased 0..2046) */
    for (size_t i = 0; i < n_random / 4; i++) {
        const uint64_t r = rng();
        buf[k++] = bitsd((r & 0x800fffffffffffffULL) | ((uint64_t)((r >> 52) % 2047) << 52));
    }
    for (size_t i = 0; i < n_random / 4; i++) buf[k++] = bitsd(rng());  /* arbitrary bit patterns */
    return k;
}
