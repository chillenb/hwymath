/* Bit-exact comparison of a C scalar callout (sin, cos, tan) against Intel's assembly routine: results,
 * whether the result was written, return codes and raised FP exception flags. Plain x86-64; runs on
 * any machine. Included by functions/{sin,cos,tan}/tests/test_rare.c, which define
 *   TRIG_NAME      "sin"
 *   TRIG_COUT_C    svml_dsin_ha_cout_rare
 *   TRIG_COUT_ASM  __svml_dsin_ha_cout_rare_internal */
#include <fenv.h>
#include <stdio.h>

#include "functions/sin/tests/inputs.h"

int TRIG_COUT_ASM(const double *px, double *pr); /* from the assembly */

#define FLAGS (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)

int
main(int argc, char **argv)
{
    size_t n_random = argc > 1 ? strtoull(argv[1], NULL, 10) : 20000000;
    double *in = malloc((n_random + 40000) * sizeof *in);
    size_t n = make_inputs(in, n_random);
    size_t bad = 0, by_code[2] = {0}, written = 0;

    for (size_t i = 0; i < n; i++) {
        double rc = 0x1.234p5, ra = 0x1.234p5; /* sentinel: the routines write only for inf/NaN */
        feclearexcept(FLAGS);
        int cc = TRIG_COUT_C(&in[i], &rc);
        int fc = fetestexcept(FLAGS);
        feclearexcept(FLAGS);
        int ca = TRIG_COUT_ASM(&in[i], &ra);
        int fa = fetestexcept(FLAGS);
        if (cc >= 0 && cc < 2) by_code[cc]++;
        written += dbits(ra) != dbits(0x1.234p5);
        if (dbits(rc) != dbits(ra) || cc != ca || fc != fa) {
            if (bad++ < 10) {
                printf("MISMATCH x=%a: C %a (code %d, flags %#x)  asm %a (code %d, flags %#x)\n", in[i], rc, cc, fc,
                       ra, ca, fa);
            }
        }
    }
    printf("test_rare (%s): %zu inputs (return codes: 0=%zu 1=%zu; result written for %zu), %zu mismatches\n",
           TRIG_NAME, n, by_code[0], by_code[1], written, bad);
    free(in);
    return bad != 0;
}
