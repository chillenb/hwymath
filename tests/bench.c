/* Throughput of exp implementations on donuts-style AVX2 hardware: 4096 random doubles in [-700, 700]
 * (stays in L1/L2), repeated; best of several runs, nanoseconds per value. */
#include <stdio.h>
#include <time.h>

#include "../svml_exp_d_ha.h"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "common.h"

#define AVX2 __attribute__((target("avx2,fma")))
AVX2 __m256d __svml_exp4_ha_l9(__m256d);
AVX2 __m256d __svml_exp4_l9(__m256d);
AVX2 __m256d _ZGVdN4v_exp(__m256d);
#ifdef HAVE_SLEEF  /* declared here: sleef.h only declares it when compiled with -mavx2 */
AVX2 __m256d Sleef_expd4_u10avx2(__m256d);
#endif

#define N 4096
static double in[N], out[N];
static volatile double sink;

typedef void (*kernel)(void);
#ifdef HAVE_SLEEF
AVX2 static void k_sleef(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, Sleef_expd4_u10avx2(_mm256_loadu_pd(in + i))); }
#endif
static void k_glibc(void) { for (int i = 0; i < N; i++) out[i] = exp(in[i]); }
AVX2 static void k_libmvec(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, _ZGVdN4v_exp(_mm256_loadu_pd(in + i))); }
AVX2 static void k_intel_ha(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, __svml_exp4_ha_l9(_mm256_loadu_pd(in + i))); }
AVX2 static void k_intel_la(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, __svml_exp4_l9(_mm256_loadu_pd(in + i))); }
static void k_port_array(void) { svml_exp_ha_array_avx2(in, out, N); }
static void k_fast_array(void) { svml_exp_fast_array_avx2(in, out, N); }
static void k_hp_array(void) { svml_exp_hp_array_avx2(in, out, N); }
AVX2 static void k_hp(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, svml_exp4_hp_avx2(_mm256_loadu_pd(in + i))); }
AVX2 static void k_fast(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, svml_exp4_fast_avx2(_mm256_loadu_pd(in + i))); }
AVX2 static void k_port(void) { for (int i = 0; i < N; i += 4) _mm256_storeu_pd(out + i, svml_exp4_ha_avx2(_mm256_loadu_pd(in + i))); }

static double
now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int
main(void)
{
    for (int i = 0; i < N; i++) in[i] = ((double)(rng() >> 11) * 0x1p-53) * 1400.0 - 700.0;
    struct { const char *name; kernel k; } ks[] = {
        {"glibc exp (scalar)", k_glibc}, {"glibc libmvec AVX2", k_libmvec}, {"Intel __svml_exp4_l9 (LA)", k_intel_la},
        {"Intel __svml_exp4_ha_l9 (HA)", k_intel_ha}, {"this port (AVX2, HA)", k_port}, {"this port, array function", k_port_array},
        {"fast variant (AVX2)", k_fast}, {"fast variant, array function", k_fast_array},
        {"hp variant (AVX2)", k_hp}, {"hp variant, array function", k_hp_array},
#ifdef HAVE_SLEEF
        {"SLEEF Sleef_expd4_u10avx2", k_sleef},
#endif
    };
    const int reps = 2000;
    for (size_t k = 0; k < sizeof ks / sizeof *ks; k++) {
        double best = 1e9;
        for (int trial = 0; trial < 7; trial++) {
            double t0 = now();
            for (int r = 0; r < reps; r++) ks[k].k();
            double dt = (now() - t0) / ((double)reps * N);
            sink = out[N - 1];
            if (dt < best) best = dt;
        }
        printf("%-30s %6.2f ns/value  %7.0f Mvalues/s\n", ks[k].name, best * 1e9, 1e-6 / best);
    }
    return 0;
}
