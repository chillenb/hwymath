/* Other acosh implementations available on this machine, as array functions, for comparison with the
 * Highway acosh (hwy_test.cc). extras_x86.c / extras_arm.c (+ extras_arm_sve.c) implement it. */
#ifndef HWYMATH_ACOSH_EXTRAS_H
#define HWYMATH_ACOSH_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*acosh_array_fn)(const double *x, double *y, size_t len);

struct extra_impl {
    const char *name;
    acosh_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
