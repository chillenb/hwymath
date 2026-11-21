/* Accuracy statistics, result hashing and timing shared by the tests of all functions. C and C++.
 *
 * Errors are measured in ulps of the correctly rounded double result against a long double reference
 * (x86-64: 64-bit mantissa, AArch64: binary128). Typical use:
 *     struct acc_stats s = ACC_STATS_INIT;
 *     for (...) acc_add(&s, x, y, ref);      // ref = sinl((long double)x) etc.
 *     acc_print(stdout, "name", &s);
 */
#ifndef HWYMATH_ACCURACY_H
#define HWYMATH_ACCURACY_H

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

/* error of y in ulps of the correctly rounded result round(ref) */
static inline double
acc_ulp_err(double y, long double ref)
{
    const double rd = (double)ref;
    long double ulp = (long double)nextafter(fabs(rd), INFINITY) - fabsl((long double)rd);
    if (fabs(rd) == DBL_MAX) ulp = 0x1p971L;
    return (double)(fabsl((long double)y - ref) / ulp);
}

struct acc_stats {
    double max_err, worst_x;  /* largest error and the (first) input where it occurs */
    size_t n, not_cr, over1;  /* counted results, not correctly rounded (> 0.5 ulp), > 1 ulp */
    size_t mismatched_special; /* NaN/inf/sign mismatches where the reference is not finite or zero */
};
#define ACC_STATS_INIT {0, 0, 0, 0, 0, 0}

/* Adds one result. Finite nonzero references count towards the ulp statistics; for non-finite or zero
 * references the result must match exactly (NaN for NaN, same infinity, same signed zero). */
static inline void
acc_add(struct acc_stats *s, double x, double y, long double ref)
{
    const double rd = (double)ref;
    if (isfinite(rd) && rd != 0.0) {
        const double e = acc_ulp_err(y, ref);
        s->n++;
        s->not_cr += e > 0.5;
        s->over1 += e > 1.0;
        if (e > s->max_err || isnan(e)) { s->max_err = isnan(e) ? INFINITY : e; s->worst_x = x; }
    } else if (isnan(rd) ? !isnan(y) : !(y == rd && signbit(y) == signbit(rd))) {
        s->mismatched_special++;
    }
}

static inline void
acc_print(FILE *f, const char *name, const struct acc_stats *s)
{
    fprintf(f, "%-26s max %.4f ulp (x = %.17g), not correctly rounded %.4f%%, > 1 ulp: %zu, special mismatches: %zu\n",
            name, s->max_err, s->worst_x, s->n ? 100.0 * s->not_cr / s->n : 0.0, s->over1, s->mismatched_special);
}

/* FNV-1a hash of the bits of n doubles (bit-identity checks across machines: make hashes) */
static inline uint64_t
acc_hash(const double *v, size_t n)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    const unsigned char *p = (const unsigned char *)v;
    for (size_t i = 0; i < n * sizeof(double); i++) h = (h ^ p[i]) * 0x100000001b3ULL;
    return h;
}

/* The same with every NaN counted as one canonical NaN, for cross-architecture comparisons of functions that
 * produce NaNs by invalid operations: the default NaN is negative on x86-64 and positive on AArch64. */
static inline uint64_t
acc_hash_nan_canonical(const double *v, size_t n)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    for (size_t k = 0; k < n; k++) {
        double d = v[k];
        if (d != d) {
            const uint64_t qnan = 0x7ff8000000000000ULL;
            memcpy(&d, &qnan, sizeof d);
        }
        const unsigned char *p = (const unsigned char *)&d;
        for (size_t i = 0; i < sizeof(double); i++) h = (h ^ p[i]) * 0x100000001b3ULL;
    }
    return h;
}

static inline double
acc_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Throughput of an array function y[i] = f(x[i]) in ns/value: best of 7 runs of reps calls over
 * the n inputs in x (keep n small enough to stay in cache, e.g. 4096). */
typedef void (*acc_array_fn)(const double *x, double *y, size_t n);
static inline double
acc_bench(acc_array_fn fn, const double *x, double *y, size_t n, int reps)
{
    double best = 1e30;
    for (int trial = 0; trial < 7; trial++) {
        const double t0 = acc_now();
        for (int r = 0; r < reps; r++) fn(x, y, n);
        const double dt = (acc_now() - t0) / ((double)reps * (double)n);
        if (dt < best) best = dt;
    }
    return best * 1e9;
}

/* Same for two-argument functions z[i] = f(x[i], y[i]) (atan2, pow). */
typedef void (*acc_array2_fn)(const double *x, const double *y, double *z, size_t n);
static inline double
acc_bench2(acc_array2_fn fn, const double *x, const double *y, double *z, size_t n, int reps)
{
    double best = 1e30;
    for (int trial = 0; trial < 7; trial++) {
        const double t0 = acc_now();
        for (int r = 0; r < reps; r++) fn(x, y, z, n);
        const double dt = (acc_now() - t0) / ((double)reps * (double)n);
        if (dt < best) best = dt;
    }
    return best * 1e9;
}

#endif
