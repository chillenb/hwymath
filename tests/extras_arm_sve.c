/* SLEEF's SVE exp, vector-length agnostic loop. Compiled with SVE enabled; only called if the CPU has SVE. */
#include <stddef.h>

#ifdef HAVE_SLEEF
#include <arm_sve.h>
#include <sleef.h>

void
sleef_sve(const double *x, double *y, size_t n)
{
    for (size_t i = 0; i < n; i += svcntd()) {
        const svbool_t pg = svwhilelt_b64((uint64_t)i, (uint64_t)n);
        svst1_f64(pg, y + i, Sleef_expdx_u10sve(svld1_f64(pg, x + i)));
    }
}
#endif
