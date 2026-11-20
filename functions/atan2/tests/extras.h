/* Other atan2 implementations available on this machine, as array functions, for comparison with the
 * Highway atan2 (hwy_test.cc): Intel's libsvml (HAVE_INTEL_SVML) and SLEEF (HAVE_SLEEF), AVX2 and AVX-512. */
#ifndef HWYMATH_ATAN2_EXTRAS_H
#define HWYMATH_ATAN2_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*atan2_array_fn)(const double *y, const double *x, double *r, size_t len);

struct extra_impl {
    const char *name;
    atan2_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
