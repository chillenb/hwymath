/* vrcp14pd (hardware or Intel SDE) vs the scalar reference rcp14_ref_bits() (Intel's published reference
 * algorithm, common/include/hwymath/rcp14_ref.h), bit for bit, on the inputs of rcp14_inputs.h. Needs AVX-512F (else exit 77). */
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>

#include "hwymath/test_util.h"
#include "hwymath/rcp14_ref.h"
#include "rcp14_inputs.h"

__attribute__((target("avx512f"))) static size_t
check(const uint64_t *in, size_t n)
{
    size_t bad = 0;
    for (size_t i = 0; i + 8 <= n; i += 8) {
        uint64_t out[8];
        _mm512_storeu_si512(out, _mm512_castpd_si512(_mm512_rcp14_pd(_mm512_loadu_pd((const double *)(in + i)))));
        for (int l = 0; l < 8; l++) {
            if (out[l] != rcp14_ref_bits(in[i + l]) && bad++ < 10) {
                printf("MISMATCH x=%016llx: vrcp14pd %016llx  reference %016llx\n", (unsigned long long)in[i + l],
                       (unsigned long long)out[l], (unsigned long long)rcp14_ref_bits(in[i + l]));
            }
        }
    }
    return bad;
}

int
main(int argc, char **argv)
{
    if (!__builtin_cpu_supports("avx512f")) {
        printf("test_rcp14: SKIPPED (this CPU has no AVX-512F)\n");
        return 77;
    }
    const size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 1000000;
    uint64_t *in = malloc((rcp14_max_inputs(n_random) + 8) * sizeof *in);
    size_t n = rcp14_inputs(in, n_random);
    while (n % 8) in[n++] = 0x3ff0000000000000ULL;
    const size_t bad = check(in, n);
    printf("test_rcp14: %zu inputs, %zu mismatches between vrcp14pd and the reference\n", n, bad);
    free(in);
    return bad != 0;
}
