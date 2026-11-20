/* Bit-exact comparison of a C/AVX-512 translation (sin, cos, tan) against Intel's assembly. Needs
 * AVX-512F/DQ at run time; exits with 77 (skipped) otherwise. Included by
 * functions/{sin,cos,tan}/tests/test_<f>8.c, which define
 *   TRIG_NAME  "sin"
 *   TRIG_C     svml_sin8_ha
 *   TRIG_ASM   __svml_sin8_ha
 * usage: test_<f>8 [N_RANDOM]   (also prints how many inputs took the large-argument path) */
#include <stdio.h>

#include "functions/sin/tests/inputs.h"

__attribute__((target("avx512f"))) __m512d TRIG_ASM(__m512d x); /* from the assembly */

__attribute__((target("avx512f,avx512dq"))) static size_t
compare(const double *in, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        double rc[8], ra[8];
        __m512d x = _mm512_loadu_pd(in + i);
        _mm512_storeu_pd(rc, TRIG_C(x));
        _mm512_storeu_pd(ra, TRIG_ASM(x));
        for (int l = 0; l < 8; l++) {
            if (dbits(rc[l]) != dbits(ra[l]) && bad++ < 10) {
                printf("MISMATCH x=%a: C %a  asm %a\n", in[i + l], rc[l], ra[l]);
            }
        }
    }
    return bad;
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512dq")) {
        printf("test_%s8: SKIPPED (this CPU has no AVX-512F/DQ)\n", TRIG_NAME);
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000;
    double *in = malloc((n_random + 40000 + 8) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    while (n % 8) in[n++] = 0.0;
    size_t big = 0;
    for (size_t i = 0; i < n; i++) big += fabs(in[i]) > TRIG_THRESHOLD && isfinite(in[i]);
    size_t bad = compare(in, n);
    printf("test_%s8: %zu inputs (%zu finite with |x| > %g: large-argument path), %zu mismatches\n", TRIG_NAME, n, big,
           TRIG_THRESHOLD, bad);
    free(in);
    return bad != 0;
}
