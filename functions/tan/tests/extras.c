/* Comparison implementations of tan for hwy_test (extras_tmpl.h). */
#include "../tan_d.h"
#define TRIG_F tan
#define TRIG_SVML_ARRAY svml_tan_ha_array
#define TRIG_NO_LIBMVEC /* not in the glibc 2.28 sysroot (libmvec has tan since glibc 2.35) */
#include "functions/sin/tests/extras_tmpl.h"
