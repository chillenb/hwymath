// Cross-platform check of the Highway log2 (static target); the program is functions/log/tests/log_family_xcheck.h.
#include <cmath>
#include <cstddef>

extern "C" {
void hwy_log2_array(const double *x, double *y, size_t len);
int hwy_log2_variant(void);
size_t hwy_log2_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log2"
#define REF(x) log2l(x)
#define HWY_ARRAY hwy_log2_array
#define VARIANT hwy_log2_variant
#define LANES hwy_log2_lanes

#include "functions/log/tests/log_family_xcheck.h"
