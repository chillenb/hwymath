/* Bit-exact comparison of svml_pow8_ha() (C intrinsics) against Intel's __svml_pow8_ha assembly.
 * Needs AVX-512F/DQ at run time; exits with 77 (skipped) otherwise. */
#include <stdio.h>

#include "../pow_d.h"
#include "inputs.h"

#define AVX512 __attribute__((target("avx512f,avx512dq")))
AVX512 __m512d __svml_pow8_ha(__m512d x, __m512d y);  /* from the assembly */

AVX512 static size_t
compare(const double *x, const double *y, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        double rc[8], ra[8];
        const __m512d vx = _mm512_loadu_pd(x + i), vy = _mm512_loadu_pd(y + i);
        _mm512_storeu_pd(rc, svml_pow8_ha(vx, vy));
        _mm512_storeu_pd(ra, __svml_pow8_ha(vx, vy));
        for (int l = 0; l < 8; l++) {
            if (dbits(rc[l]) != dbits(ra[l]) && bad++ < 10) {
                printf("MISMATCH x=%a y=%a: C %a  asm %a\n", x[i + l], y[i + l], rc[l], ra[l]);
            }
        }
    }
    return bad;
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512dq")) {
        printf("test_pow8: SKIPPED (this CPU has no AVX-512F/DQ)\n");
        return 77;
    }
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000;
    double *x = malloc((n_random + POW_INPUTS_FIXED + 8) * sizeof *x);
    double *y = malloc((n_random + POW_INPUTS_FIXED + 8) * sizeof *y);
    size_t n = make_inputs(x, y, n_random);
    while (n % 8) { x[n] = 1.0; y[n] = 1.0; n++; }
    size_t bad = compare(x, y, n);
    printf("test_pow8: %zu inputs, %zu mismatches\n", n, bad);
    free(x);
    free(y);
    return bad != 0;
}
