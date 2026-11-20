/* Checks the model of SVML cbrt's reciprocal step (rcp14_model.h) against the instructions it replaces,
 * vrndscalepd(vrcp14pd(m), 0x48): exhaustively over all 2^16 16-bit mantissa prefixes of m in [1, 2)
 * (both signs; low 36 bits zero, all ones and random), then for random m in +-[1, 2).
 * Also verifies that vrcp14pd itself depends only on those top 16 bits (unrounded output equal for all
 * low-bit patterns of a prefix), except that rcp14(+-1) = +-1 exactly. Needs AVX-512F at run time; exits with 77 (skipped) otherwise. */
#include <stdio.h>
#include <stdlib.h>

#include "hwymath/test_util.h"
#include "../rcp14_model.h"
#include <immintrin.h>

__attribute__((target("avx512f"))) static void
rcp_rounded(const double *x, double *rcp, double *rounded)
{
    const __m512d r = _mm512_rcp14_pd(_mm512_loadu_pd(x));
    _mm512_storeu_pd(rcp, r);
    _mm512_storeu_pd(rounded, _mm512_roundscale_round_pd(r, 0x48, _MM_FROUND_NO_EXC));
}

static size_t
check(const double x[8], size_t *bad_prefix)
{
    double rcp[8], rr[8];
    size_t bad = 0;
    rcp_rounded(x, rcp, rr);
    for (int l = 0; l < 8; l++) {
        const int k = cbrt_rcp14_k(x[l]);
        const double model = (x[l] < 0 ? -k : k) / 16.0;
        if (dbits(model) != dbits(rr[l]) && bad++ < 10) {
            printf("MISMATCH m=%a: rcp14 %a rounded %a, model %a\n", x[l], rcp[l], rr[l], model);
        }
        /* rcp14 depends on the top 16 mantissa bits only: compare with the prefix's rcp14 (lane 1, low bits
         * all ones); m = +-1 (exact power of two) gives exactly +-1 */
        const double ref = (x[l] == 1.0 || x[l] == -1.0) ? 1.0 : (rcp[1] < 0 ? -rcp[1] : rcp[1]);
        if (bad_prefix && dbits(rcp[l] < 0 ? -rcp[l] : rcp[l]) != dbits(ref)) {
            if ((*bad_prefix)++ < 10) printf("rcp14 depends on low bits: m=%a\n", x[l]);
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
    const size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 4000000;
    size_t bad = 0, bad_prefix = 0, n = 0;
    for (uint64_t p = 0; p < 65536; p++) {
        const uint64_t hi = 0x3ff0000000000000ULL | (p << 36);
        double x[8];
        x[0] = bitsd(hi);
        x[1] = bitsd(hi | 0xfffffffffULL);
        x[2] = bitsd(hi | (rng() & 0xfffffffffULL));
        x[3] = bitsd(hi | 1);
        for (int l = 0; l < 4; l++) x[4 + l] = -x[l];
        bad += check(x, &bad_prefix);
        n += 8;
    }
    for (size_t i = 0; i < n_random; i += 8) {
        double x[8];
        for (int l = 0; l < 8; l++) x[l] = bitsd((rng() & 0x800fffffffffffffULL) | 0x3ff0000000000000ULL);
        bad += check(x, NULL);
        n += 8;
    }
    printf("test_rcp14: %zu inputs (all 2^16 prefixes), %zu mismatches of the model, "
           "%zu cases where rcp14 depends on more than 16 mantissa bits\n", n, bad, bad_prefix);
    return bad != 0 || bad_prefix != 0;
}
