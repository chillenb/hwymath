/* Bit-exact comparison of svml_atan28_ha() (C intrinsics) against Intel's __svml_atan28_ha assembly
 * (y in zmm0, x in zmm1). Needs AVX-512F at run time; exits with 77 (skipped) otherwise. */
#include <stdio.h>

#include "../atan2_d.h"
#include "inputs.h"

__attribute__((target("avx512f"))) __m512d __svml_atan28_ha(__m512d y, __m512d x);  /* from the assembly */

__attribute__((target("avx512f"))) static size_t
compare(const double *ys, const double *xs, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        double rc[8], ra[8];
        const __m512d y = _mm512_loadu_pd(ys + i), x = _mm512_loadu_pd(xs + i);
        _mm512_storeu_pd(rc, svml_atan28_ha(y, x));
        _mm512_storeu_pd(ra, __svml_atan28_ha(y, x));
        for (int l = 0; l < 8; l++) {
            if (dbits(rc[l]) != dbits(ra[l]) && bad++ < 10) {
                printf("MISMATCH y=%a x=%a: C %a  asm %a\n", ys[i + l], xs[i + l], rc[l], ra[l]);
            }
        }
    }
    return bad;
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx512f")) {
        printf("test_atan28: SKIPPED (this CPU has no AVX-512F)\n");
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000;
    double *ys = malloc((n_random + INPUTS_EXTRA + 8) * sizeof *ys);
    double *xs = malloc((n_random + INPUTS_EXTRA + 8) * sizeof *xs);
    size_t n = make_inputs(ys, xs, n_random);
    while (n % 8) { ys[n] = xs[n] = 1.0; n++; }
    size_t bad = compare(ys, xs, n);
    printf("test_atan28: %zu input pairs, %zu mismatches\n", n, bad);
    free(ys);
    free(xs);
    return bad != 0;
}
