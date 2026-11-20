/* Compare the AVX2 port with Intel's AVX2 SVML exp (libsvml.so from the Intel runtime: __svml_exp4_ha_l9,
 * __svml_exp4_l9), glibc's exp and libmvec's AVX2 exp, against a long double reference (expl).
 * Reports bit-identity with Intel's HA version and each implementation's error in ulps. */
#include <float.h>
#include <stdio.h>

#include "../svml_exp_d_ha.h"
#include "common.h"

#define AVX2 __attribute__((target("avx2,fma")))
AVX2 __m256d __svml_exp4_ha_l9(__m256d);  /* Intel libsvml, AVX2 high accuracy */
AVX2 __m256d __svml_exp4_l9(__m256d);     /* Intel libsvml, AVX2 default (lower) accuracy */
AVX2 __m256d _ZGVdN4v_exp(__m256d);       /* glibc libmvec, AVX2 */
#ifdef HAVE_SLEEF  /* declared here: sleef.h only declares it when compiled with -mavx2 */
AVX2 __m256d Sleef_expd4_u10avx2(__m256d);
#endif

#ifdef HAVE_SLEEF
enum { MINE, FAST, HP, INTEL_HA, INTEL_LA, GLIBC, LIBMVEC, SLEEF, NIMPL };
#else
enum { MINE, FAST, HP, INTEL_HA, INTEL_LA, GLIBC, LIBMVEC, NIMPL };
#endif
static const char *names[NIMPL] = {"this port (AVX2)", "fast variant (AVX2)", "hp variant (AVX2)", "Intel __svml_exp4_ha_l9", "Intel __svml_exp4_l9",
                                   "glibc exp", "glibc libmvec AVX2",
#ifdef HAVE_SLEEF
                                   "SLEEF expd4_u10avx2",
#endif
};

AVX2 static void
eval4(const double *x, double out[NIMPL][4])
{
    __m256d v = _mm256_loadu_pd(x);
    _mm256_storeu_pd(out[MINE], svml_exp4_ha_avx2(v));
    _mm256_storeu_pd(out[FAST], svml_exp4_fast_avx2(v));
    _mm256_storeu_pd(out[HP], svml_exp4_hp_avx2(v));
    _mm256_storeu_pd(out[INTEL_HA], __svml_exp4_ha_l9(v));
    _mm256_storeu_pd(out[INTEL_LA], __svml_exp4_l9(v));
    _mm256_storeu_pd(out[LIBMVEC], _ZGVdN4v_exp(v));
#ifdef HAVE_SLEEF
    _mm256_storeu_pd(out[SLEEF], Sleef_expd4_u10avx2(v));
#endif
    for (int l = 0; l < 4; l++) out[GLIBC][l] = exp(x[l]);
}

/* error of y in ulps of the correctly rounded result; ref is the long double reference */
static double
ulp_err(double y, long double ref)
{
    double rd = (double)ref;
    long double ulp = (long double)nextafter(fabs(rd), INFINITY) - fabsl((long double)rd);
    if (fabs(rd) == DBL_MAX) ulp = 0x1p971L;
    return (double)(fabsl((long double)y - ref) / ulp);
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx2") || !__builtin_cpu_supports("fma")) {
        printf("needs AVX2 and FMA\n");
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *in = malloc((n_random + 4096 + 4) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    while (n % 4) in[n++] = 0.0;

    size_t same_ha = 0, diff_ha = 0, n_stat = 0, n_special = 0, special_diff = 0;
    double max_err[NIMPL] = {0}, worst_x[NIMPL] = {0};
    size_t over_half[NIMPL] = {0}, over_one[NIMPL] = {0};
    for (size_t i = 0; i < n; i += 4) {
        double out[NIMPL][4];
        eval4(in + i, out);
        for (int l = 0; l < 4; l++) {
            double x = in[i + l];
            int same = dbits(out[MINE][l]) == dbits(out[INTEL_HA][l]);
            same_ha += same;
            if (!same && diff_ha++ < 8) {
                printf("differs from Intel HA: x=%a (%.17g): port %a, Intel %a\n", x, x, out[MINE][l], out[INTEL_HA][l]);
            }
            /* special inputs/results: the fast variant must agree exactly with the HA port */
            double h = out[MINE][l];
            if (!isfinite(x) || !isfinite(h) || h == 0.0 || fabs(h) < 0x1p-1022) {
                n_special++;
                if ((dbits(out[FAST][l]) != dbits(h) || dbits(out[HP][l]) != dbits(h)) && special_diff++ < 8) {
                    printf("special input differs: x=%a: fast %a, hp %a, HA port %a\n", x, out[FAST][l], out[HP][l], h);
                }
            }
            /* error statistics where the result is a finite, nonzero double */
            long double ref = expl((long double)x);
            if (!isfinite(x) || !isfinite((double)ref) || (double)ref == 0.0) continue;
            n_stat++;
            for (int k = 0; k < NIMPL; k++) {
                double e = ulp_err(out[k][l], ref);
                if (e > max_err[k]) { max_err[k] = e; worst_x[k] = x; }
                over_half[k] += e > 0.5;
                over_one[k] += e > 1.0;
            }
        }
    }
    printf("\n%zu inputs: port vs Intel __svml_exp4_ha_l9 bit-identical for %zu, different for %zu\n\n",
           n, same_ha, diff_ha);
    printf("special inputs/results (NaN, inf, 0, denormal): %zu, fast/hp variants differ from the HA port on %zu\n\n",
           n_special, special_diff);
    printf("error vs expl() reference over %zu finite results:\n", n_stat);
    printf("%-26s %10s %14s %12s   %s\n", "implementation", "max ulp", "not corr.rnd.", "> 1 ulp", "worst x");
    for (int k = 0; k < NIMPL; k++) {
        printf("%-26s %10.4f %13.5f%% %12zu   %.17g\n", names[k], max_err[k],
               100.0 * over_half[k] / n_stat, over_one[k], worst_x[k]);
    }
    free(in);
    return 0;
}
