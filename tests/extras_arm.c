/* AArch64 comparison implementations: SLEEF Advanced SIMD and, if the CPU has SVE, SLEEF SVE (extras_arm_sve.c). */
#include <string.h>
#include <sys/auxv.h>

#include "extras.h"

#ifndef HWCAP_SVE
#define HWCAP_SVE (1 << 22)
#endif

#ifdef HAVE_SLEEF
#include <arm_neon.h>
#include <sleef.h>

static void
sleef_advsimd(const double *x, double *y, size_t n)
{
    size_t i = 0;
    for (; i + 2 <= n; i += 2) vst1q_f64(y + i, Sleef_expd2_u10advsimd(vld1q_f64(x + i)));
    for (; i < n; i++) {
        double b[2] = {x[i], 0}, o[2];
        vst1q_f64(o, Sleef_expd2_u10advsimd(vld1q_f64(b)));
        y[i] = o[0];
    }
}

void sleef_sve(const double *x, double *y, size_t n);  /* extras_arm_sve.c */
#endif

size_t
extra_impls(struct extra_impl *out, size_t max)
{
    size_t k = 0;
#ifdef HAVE_SLEEF
    if (k < max) { out[k].name = "SLEEF u10 AdvSIMD"; out[k].fn = sleef_advsimd; k++; }
    if ((getauxval(AT_HWCAP) & HWCAP_SVE) && k < max) { out[k].name = "SLEEF u10 SVE"; out[k].fn = sleef_sve; k++; }
#else
    (void)out; (void)max;
#endif
    return k;
}
