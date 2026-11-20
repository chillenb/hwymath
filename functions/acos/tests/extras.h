/* Other acos implementations available on this machine, as array functions, for comparison with the
 * Highway acos (hwy_test.cc): Intel's libsvml (HAVE_INTEL_SVML) and SLEEF (HAVE_SLEEF), AVX2 and AVX-512. */
#ifndef HWYMATH_ACOS_EXTRAS_H
#define HWYMATH_ACOS_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*acos_array_fn)(const double *x, double *y, size_t len);

struct extra_impl {
    const char *name;
    acos_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
