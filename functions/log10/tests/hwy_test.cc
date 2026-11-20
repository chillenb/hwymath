// Test and benchmark of the Highway log10 (hwy/hwy_log10.cc); the test itself is functions/log/tests/log_family_test.h.
#include <cmath>
#include <cstddef>

extern "C" {
#include "../log10_d.h"
#include "functions/log/tests/extras.h"
void hwy_log10_array(const double *x, double *y, size_t len);
void hwy_contrib_log10_array(const double *x, double *y, size_t len);
int hwy_log10_variant(void);
size_t hwy_log10_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log10"
#define REF(x) log10l(x)
#define LIBM(x) std::log10(x)
#define HWY_ARRAY hwy_log10_array
#define CONTRIB_ARRAY hwy_contrib_log10_array
#define VARIANT hwy_log10_variant
#define LANES hwy_log10_lanes
#define CALLOUT svml_dlog10_ha_cout_rare
#define IS_SPECIAL(x) (!((x) > 0.0 && (x) < HUGE_VAL))
#define SVML_ARRAY svml_log10_ha_array
#define SVML_ISA_OK() (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
#define BENCH_LO 0.001
#define BENCH_HI 1000.0

#include "functions/log/tests/log_family_test.h"
