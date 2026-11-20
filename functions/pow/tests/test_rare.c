/* Bit-exact comparison of svml_dpow_ha_cout_rare() (C) against Intel's scalar assembly routine
 * __svml_dpow_ha_cout_rare_internal: results, return codes and raised FP exception flags, for every input
 * pair (the callout handles any (x, y), not only the lanes the vector path flags). Plain x86-64. */
#include <fenv.h>
#include <stdio.h>

#include "../pow_d.h"
#include "inputs.h"

int __svml_dpow_ha_cout_rare_internal(const double *px, const double *py, double *pr);  /* from the assembly */

#define FLAGS (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *x = malloc((n_random + POW_INPUTS_FIXED) * sizeof *x);
    double *y = malloc((n_random + POW_INPUTS_FIXED) * sizeof *y);
    size_t n = make_inputs(x, y, n_random);
    size_t bad = 0, by_code[2] = {0}, flagged[5] = {0};

    for (size_t i = 0; i < n; i++) {
        double rc = 0, ra = 0;
        feclearexcept(FLAGS);
        int cc = svml_dpow_ha_cout_rare(&x[i], &y[i], &rc);
        int fc = fetestexcept(FLAGS);
        feclearexcept(FLAGS);
        int ca = __svml_dpow_ha_cout_rare_internal(&x[i], &y[i], &ra);
        int fa = fetestexcept(FLAGS);
        if (cc >= 0 && cc < 2) by_code[cc]++;
        flagged[0] += !!(fa & FE_INEXACT); flagged[1] += !!(fa & FE_UNDERFLOW); flagged[2] += !!(fa & FE_OVERFLOW);
        flagged[3] += !!(fa & FE_INVALID); flagged[4] += !!(fa & FE_DIVBYZERO);
        if (dbits(rc) != dbits(ra) || cc != ca || fc != fa) {
            if (bad++ < 10) {
                printf("MISMATCH x=%a y=%a: C %a (code %d, flags %#x)  asm %a (code %d, flags %#x)\n",
                       x[i], y[i], rc, cc, fc, ra, ca, fa);
            }
        }
    }
    printf("test_rare: %zu inputs (return codes: 0=%zu 1=%zu; flags raised: inexact %zu, underflow %zu, "
           "overflow %zu, invalid %zu, divbyzero %zu), %zu mismatches\n",
           n, by_code[0], by_code[1], flagged[0], flagged[1], flagged[2], flagged[3], flagged[4], bad);
    free(x);
    free(y);
    return bad != 0;
}
