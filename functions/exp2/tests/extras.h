/* Other exp2 implementations available on this machine, as array functions, for comparison with the
 * Highway exp2 (hwy_test.cc). Implemented in extras.c (x86-64 only; empty elsewhere). */
#ifndef HWYMATH_EXP2_EXTRAS_H
#define HWYMATH_EXP2_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*array_fn)(const double *x, double *y, size_t len);

struct extra_impl {
    const char *name;
    array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
