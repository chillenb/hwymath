// Cross-platform check of the Highway log (static target); the program is functions/log/tests/log_family_xcheck.h.
#include <cmath>
#include <cstddef>

extern "C" {
void hwy_log_array(const double *x, double *y, size_t len);
int hwy_log_variant(void);
size_t hwy_log_lanes(void);
}
#include "inputs.h"

#define FN_NAME "log"
#define REF(x) logl(x)
#define HWY_ARRAY hwy_log_array
#define VARIANT hwy_log_variant
#define LANES hwy_log_lanes

#include "functions/log/tests/log_family_xcheck.h"
