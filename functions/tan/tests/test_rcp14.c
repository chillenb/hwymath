/* vrcp14pd (hardware or Intel SDE) vs the scalar reference rcp14_ref_bits() (Intel's published reference
 * algorithm, tests/rcp14_ref.h), bit for bit: every 16-bit leading fraction pattern (the instruction only
 * reads 16 fraction bits) with random low bits, at the exponents with special handling and some
 * others, both signs; denormals; special values; random bit patterns. Needs AVX-512F (else exit 77). */
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>

#include "hwymath/test_util.h"
#include "rcp14_ref.h"

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
    static const uint64_t exps[] = {0x3ff, 0x3fe, 0x400, 0x001, 0x002, 0x7fc, 0x7fd, 0x7fe, 0x123, 0x5a5};
    const size_t nexp = sizeof exps / sizeof *exps;
    uint64_t *in = malloc((65536 * (nexp + 1) * 2 + n_random + 64) * sizeof *in);
    size_t n = 0;
    for (size_t k = 0; k < nexp; k++) {
        for (uint64_t f = 0; f < 65536; f++) {
            const uint64_t low = (k & 1) ? rng() & 0xfffffffffULL : 0;  /* the 36 ignored bits */
            in[n++] = exps[k] << 52 | f << 36 | low;
            in[n++] = 0x8000000000000000ULL | exps[k] << 52 | f << 36 | (rng() & 0xfffffffffULL);
        }
    }
    for (uint64_t f = 0; f < 65536; f++) {  /* denormals with every leading 16-bit pattern after normalization */
        const int sh = (int)(f % 52);
        in[n++] = ((0x0010000000000000ULL | f << 36 | (rng() & 0xfffffffffULL)) >> (sh + 1));
        in[n++] = 0x8000000000000000ULL | (rng() >> (12 + sh));
    }
    const uint64_t specials[] = {0, 0x8000000000000000ULL, 0x7ff0000000000000ULL, 0xfff0000000000000ULL,
                                 0x7ff8000000000000ULL, 0x7ff0000000000001ULL, 0xfff4000000000123ULL,
                                 0x0004000000000000ULL, 0x0004000000000001ULL, 0x8004000000000000ULL,
                                 0x0000000000000001ULL, 0x000fffffffffffffULL, 0x7fefffffffffffffULL,
                                 0x3ff0000000000000ULL, 0x7fe0000000000000ULL, 0x7fd0000000000000ULL};
    for (size_t i = 0; i < sizeof specials / sizeof *specials; i++) in[n++] = specials[i];
    for (size_t i = 0; i < n_random; i++) in[n++] = rng();
    while (n % 8) in[n++] = 0x3ff0000000000000ULL;
    const size_t bad = check(in, n);
    printf("test_rcp14: %zu inputs, %zu mismatches between vrcp14pd and the reference\n", n, bad);
    free(in);
    return bad != 0;
}
