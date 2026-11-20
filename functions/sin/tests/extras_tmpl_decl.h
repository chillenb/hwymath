/* Interface of the comparison implementations (extras_tmpl.h) used by hwy_test_tmpl.h. C and C++. */
#ifndef HWYMATH_TRIG_EXTRAS_H
#define HWYMATH_TRIG_EXTRAS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*trig_array_fn)(const double *x, double *y, size_t len);

struct extra_impl {
    const char *name;
    trig_array_fn fn;
};

/* fills out[] with the implementations usable on this CPU (and built in); returns how many */
size_t extra_impls(struct extra_impl *out, size_t max);

#ifdef __cplusplus
}
#endif

#endif
