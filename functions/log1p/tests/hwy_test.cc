// Test and benchmark of the Highway log1p (hwy/hwy_log1p.cc); the test itself is functions/log/tests/log_family_test.h.
#include <cmath>
#include <cstddef>

extern "C" {
#include "../log1p_d.h"
#include "functions/log/tests/extras.h"
void hwy_log1p_array(const double *x, double *y, size_t len);
void hwy_contrib_log1p_array(const double *x, double *y, size_t len);
int hwy_log1p_variant(void);
size_t hwy_log1p_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log1p"
#define REF(x) log1pl(x)
#define LIBM(x) std::log1p(x)
#define HWY_ARRAY hwy_log1p_array
#define CONTRIB_ARRAY hwy_contrib_log1p_array
#define VARIANT hwy_log1p_variant
#define LANES hwy_log1p_lanes
#define CALLOUT svml_dlog1p_ha_cout_rare
#define IS_SPECIAL(x) (!((x) + 1.0 > 0.0 && (x) + 1.0 < HUGE_VAL))
#define SVML_ARRAY svml_log1p_ha_array
#define SVML_ISA_OK() __builtin_cpu_supports("avx512f")
#define BENCH_LO -0.9
#define BENCH_HI 10.0

#include "functions/log/tests/log_family_test.h"
