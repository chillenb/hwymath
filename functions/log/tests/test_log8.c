/* Bit-exact comparison of svml_log8_ha() (C intrinsics) against Intel's __svml_log8_ha assembly.
 * Needs AVX-512F/DQ at run time; exits with 77 (skipped) otherwise. */
#include <stdio.h>

#include "../log_d.h"
#include "inputs.h"

__attribute__((target("avx512f"))) __m512d __svml_log8_ha(__m512d x);  /* from the assembly */

__attribute__((target("avx512f,avx512dq"))) static size_t
compare(const double *in, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        double rc[8], ra[8];
        __m512d x = _mm512_loadu_pd(in + i);
        _mm512_storeu_pd(rc, svml_log8_ha(x));
        _mm512_storeu_pd(ra, __svml_log8_ha(x));
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
        printf("test_log8: SKIPPED (this CPU has no AVX-512F/DQ)\n");
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000;
    double *in = malloc((n_random + 40000 + 8) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    while (n % 8) in[n++] = 1.0;
    size_t bad = compare(in, n);
    printf("test_log8: %zu inputs, %zu mismatches\n", n, bad);
    free(in);
    return bad != 0;
}
