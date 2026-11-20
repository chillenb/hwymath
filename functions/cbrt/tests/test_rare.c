/* Bit-exact comparison of svml_dcbrt_ha_cout_rare() (C) against Intel's scalar assembly routine
 * __svml_dcbrt_ha_cout_rare_internal: results, return codes and raised FP exception flags.
 * Plain x86-64; runs on any machine. */
#include <fenv.h>
#include <stdio.h>

#include "../cbrt_d.h"
#include "inputs.h"

int __svml_dcbrt_ha_cout_rare_internal(const double *px, double *pr);  /* from the assembly */

#define FLAGS (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *in = malloc((n_random + CBRT_FIXED_INPUTS) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    size_t bad = 0, nonzero_code = 0;

    for (size_t i = 0; i < n; i++) {
        double rc = 0, ra = 0;
        feclearexcept(FLAGS);
        int cc = svml_dcbrt_ha_cout_rare(&in[i], &rc);
        int fc = fetestexcept(FLAGS);
        feclearexcept(FLAGS);
        int ca = __svml_dcbrt_ha_cout_rare_internal(&in[i], &ra);
        int fa = fetestexcept(FLAGS);
        nonzero_code += ca != 0;
        if (dbits(rc) != dbits(ra) || cc != ca || fc != fa) {
            if (bad++ < 10) {
                printf("MISMATCH x=%a: C %a (code %d, flags %#x)  asm %a (code %d, flags %#x)\n",
                       in[i], rc, cc, fc, ra, ca, fa);
            }
        }
    }
    printf("test_rare: %zu inputs (non-zero return codes: %zu), %zu mismatches\n", n, nonzero_code, bad);
    free(in);
    return bad != 0;
}
