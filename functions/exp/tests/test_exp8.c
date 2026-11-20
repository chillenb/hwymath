/* Bit-exact comparison of svml_exp8_ha() (C intrinsics) against Intel's __svml_exp8_ha assembly.
 * Needs AVX-512F at run time; exits with 77 (skipped) otherwise. */
#include <stdio.h>

#include "../exp_d.h"
#include "inputs.h"

__attribute__((target("avx512f"))) __m512d __svml_exp8_ha(__m512d x);  /* from the assembly */

__attribute__((target("avx512f"))) static size_t
compare(const double *in, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        double rc[8], ra[8];
        __m512d x = _mm512_loadu_pd(in + i);
        double r2[8];
        _mm512_storeu_pd(rc, svml_exp8_ha(x));
        _mm512_storeu_pd(ra, __svml_exp8_ha(x));
        _mm256_storeu_pd(r2, svml_exp4_ha_avx2(_mm256_loadu_pd(in + i)));      /* AVX2 port too */
        _mm256_storeu_pd(r2 + 4, svml_exp4_ha_avx2(_mm256_loadu_pd(in + i + 4)));
        for (int l = 0; l < 8; l++) {
            if (dbits(rc[l]) != dbits(ra[l]) && bad++ < 10) {
                printf("MISMATCH x=%a: C %a  asm %a\n", in[i + l], rc[l], ra[l]);
            }
            if (dbits(r2[l]) != dbits(ra[l]) && bad++ < 10) {
                printf("MISMATCH (AVX2 port) x=%a: AVX2 %a  asm %a\n", in[i + l], r2[l], ra[l]);
            }
        }
    }
    return bad;
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx512f")) {
        printf("test_exp8: SKIPPED (this CPU has no AVX-512F)\n");
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000;
    double *in = malloc((n_random + 4096 + 8) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    while (n % 8) in[n++] = 0.0;
    size_t bad = compare(in, n);
    printf("test_exp8: %zu inputs, %zu mismatches\n", n, bad);
    free(in);
    return bad != 0;
}
