/* Other pow implementations available on this machine, as array functions z[i] = pow(x[i], y[i]), for
 * comparison with the Highway pow (hwy_test.cc). extras_x86.c / extras_arm.c implement it. */
#ifndef HWYMATH_POW_EXTRAS_H
#define HWYMATH_POW_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*pow_array_fn)(const double *x, const double *y, double *z, size_t len);

struct extra_impl {
    const char *name;
    pow_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
