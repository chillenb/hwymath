/* AArch64 comparison implementations: SLEEF Advanced SIMD. */
#include <string.h>

#include "extras.h"

#ifdef HAVE_SLEEF
#include <arm_neon.h>
#include <sleef.h>

static void
sleef_advsimd(const double *x, const double *y, double *z, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) vst1q_f64(z + i, Sleef_powd2_u10advsimd(vld1q_f64(x + i), vld1q_f64(y + i)));
    for (; i < n; i++) {
        double a[2] = {x[i], 1}, b[2] = {y[i], 1}, o[2];
        vst1q_f64(o, Sleef_powd2_u10advsimd(vld1q_f64(a), vld1q_f64(b)));
        z[i] = o[0];
    }
}
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#ifdef HAVE_SLEEF
    if (k < max) { out[k].name = "SLEEF u10 AdvSIMD"; out[k].fn = sleef_advsimd; k++; }
#else
    (void)out; (void)max;
#endif
    return k;
}
