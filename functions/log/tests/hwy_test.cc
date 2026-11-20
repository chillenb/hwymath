// Test and benchmark of the Highway log (hwy/hwy_log.cc); the test itself is functions/log/tests/log_family_test.h.
#include <cmath>
#include <cstddef>

extern "C" {
#include "../log_d.h"
#include "extras.h"
void hwy_log_array(const double *x, double *y, size_t len);
void hwy_contrib_log_array(const double *x, double *y, size_t len);
int hwy_log_variant(void);
size_t hwy_log_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log"
#define REF(x) logl(x)
#define LIBM(x) std::log(x)
#define HWY_ARRAY hwy_log_array
#define CONTRIB_ARRAY hwy_contrib_log_array
#define VARIANT hwy_log_variant
#define LANES hwy_log_lanes
#define CALLOUT svml_dlog_ha_cout_rare
#define IS_SPECIAL(x) (!((x) > 0.0 && (x) < HUGE_VAL))
#define SVML_ARRAY svml_log_ha_array
#define SVML_ISA_OK() (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq"))
#define BENCH_LO 0.001
#define BENCH_HI 1000.0

#include "functions/log/tests/log_family_test.h"
