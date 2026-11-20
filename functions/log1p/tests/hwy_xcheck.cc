// Cross-platform check of the Highway log1p (static target); the program is functions/log/tests/log_family_xcheck.h.
#include <cmath>
#include <cstddef>

extern "C" {
void hwy_log1p_array(const double *x, double *y, size_t len);
int hwy_log1p_variant(void);
size_t hwy_log1p_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log1p"
#define REF(x) log1pl(x)
#define HWY_ARRAY hwy_log1p_array
#define VARIANT hwy_log1p_variant
#define LANES hwy_log1p_lanes

#include "functions/log/tests/log_family_xcheck.h"
