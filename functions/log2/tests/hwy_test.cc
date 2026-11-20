// Test and benchmark of the Highway log2 (hwy/hwy_log2.cc); the test itself is functions/log/tests/log_family_test.h.
#include <cmath>
#include <cstddef>

extern "C" {
#include "../log2_d.h"
#include "functions/log/tests/extras.h"
void hwy_log2_array(const double *x, double *y, size_t len);
void hwy_contrib_log2_array(const double *x, double *y, size_t len);
int hwy_log2_variant(void);
size_t hwy_log2_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log2"
#define REF(x) log2l(x)
#define LIBM(x) std::log2(x)
#define HWY_ARRAY hwy_log2_array
#define CONTRIB_ARRAY hwy_contrib_log2_array
#define VARIANT hwy_log2_variant
#define LANES hwy_log2_lanes
#define CALLOUT svml_dlog2_ha_cout_rare
#define IS_SPECIAL(x) (!((x) > 0.0 && (x) < HUGE_VAL))
#define SVML_ARRAY svml_log2_ha_array
#define SVML_ISA_OK() (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
#define BENCH_LO 0.001
#define BENCH_HI 1000.0

#include "functions/log/tests/log_family_test.h"
