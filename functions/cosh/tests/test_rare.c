/* Bit-exact comparison of svml_dcosh_ha_cout_rare() (C) against Intel's scalar assembly routine
 * __svml_dcosh_ha_cout_rare_internal: results, return codes and raised FP exception flags.
 * Plain x86-64; runs on any machine. */
#include <fenv.h>
#include <stdio.h>

#include "../cosh_d.h"
#include "inputs.h"

int __svml_dcosh_ha_cout_rare_internal(const double *px, double *pr);  /* from the assembly */

#define FLAGS (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *in = malloc((n_random + 8192) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    size_t bad = 0, by_code[4] = {0};

    for (size_t i = 0; i < n; i++) {
        double rc = 0, ra = 0;
        feclearexcept(FLAGS);
        int cc = svml_dcosh_ha_cout_rare(&in[i], &rc);
        int fc = fetestexcept(FLAGS);
        feclearexcept(FLAGS);
        int ca = __svml_dcosh_ha_cout_rare_internal(&in[i], &ra);
        int fa = fetestexcept(FLAGS);
        if (cc >= 0 && cc < 4) by_code[cc]++;
        if (dbits(rc) != dbits(ra) || cc != ca || fc != fa) {
            if (bad++ < 10) {
                printf("MISMATCH x=%a: C %a (code %d, flags %#x)  asm %a (code %d, flags %#x)\n",
                       in[i], rc, cc, fc, ra, ca, fa);
            }
        }
    }
    printf("test_rare: %zu inputs (return codes: 0=%zu 3=%zu), %zu mismatches\n", n, by_code[0], by_code[3], bad);
    free(in);
    return bad != 0;
}
