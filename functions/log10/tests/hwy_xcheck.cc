// Cross-platform check of the Highway log10 (static target); the program is functions/log/tests/log_family_xcheck.h.
#include <cmath>
#include <cstddef>

extern "C" {
void hwy_log10_array(const double *x, double *y, size_t len);
int hwy_log10_variant(void);
size_t hwy_log10_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log10"
#define REF(x) log10l(x)
#define HWY_ARRAY hwy_log10_array
#define VARIANT hwy_log10_variant
#define LANES hwy_log10_lanes

#include "functions/log/tests/log_family_xcheck.h"
