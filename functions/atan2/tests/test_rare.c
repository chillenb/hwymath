/* Bit-exact comparison of svml_datan2_ha_cout_rare() (C) against Intel's scalar assembly routine
 * __svml_datan2_ha_cout_rare_internal: results, return codes and raised FP exception flags.
 * Plain x86-64; runs on any machine. */
#include <fenv.h>
#include <stdio.h>

#include "../atan2_d.h"
#include "inputs.h"

int __svml_datan2_ha_cout_rare_internal(const double *py, const double *px, double *pr);  /* from the assembly */

#define FLAGS (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *ys = malloc((n_random + INPUTS_EXTRA) * sizeof *ys);
    double *xs = malloc((n_random + INPUTS_EXTRA) * sizeof *xs);
    size_t n = make_inputs(ys, xs, n_random);
    size_t bad = 0, nonzero_code = 0;

    for (size_t i = 0; i < n; i++) {
        double rc = 0, ra = 0;
        feclearexcept(FLAGS);
        int cc = svml_datan2_ha_cout_rare(&ys[i], &xs[i], &rc);
        int fc = fetestexcept(FLAGS);
        feclearexcept(FLAGS);
        int ca = __svml_datan2_ha_cout_rare_internal(&ys[i], &xs[i], &ra);
        int fa = fetestexcept(FLAGS);
        nonzero_code += cc != 0;
        if (dbits(rc) != dbits(ra) || cc != ca || fc != fa) {
            if (bad++ < 10) {
                printf("MISMATCH y=%a x=%a: C %a (code %d, flags %#x)  asm %a (code %d, flags %#x)\n",
                       ys[i], xs[i], rc, cc, fc, ra, ca, fa);
            }
        }
    }
    printf("test_rare: %zu input pairs (nonzero return codes: %zu), %zu mismatches\n", n, nonzero_code, bad);
    free(ys);
    free(xs);
    return bad != 0;
}
