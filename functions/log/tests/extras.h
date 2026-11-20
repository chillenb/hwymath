/* Other implementations of a log-family function available on this machine, as array functions, for
 * comparison with the Highway version (log_family_test.h). Shared by functions/{log,log2,log10,log1p}:
 * extras_x86.c / extras_arm.c (+ extras_arm_sve.c) are compiled once per function with -DFN=<function>. */
#ifndef HWYMATH_LOG_EXTRAS_H
#define HWYMATH_LOG_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*log_array_fn)(const double *x, double *y, size_t len);

struct extra_impl {
    const char *name;
    log_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
