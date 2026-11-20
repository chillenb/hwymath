// Highway sin built for one static target: results to files, accuracy and special values (hwy_xcheck_tmpl.h).
#include <cmath>
#define TRIG_NAME "sin"
#define TRIG_HWY_ARRAY hwy_sin_array
#define TRIG_HWY_VARIANT hwy_sin_variant
#define TRIG_HWY_LANES hwy_sin_lanes
#define TRIG_REFL sinl
#include "functions/sin/tests/hwy_xcheck_tmpl.h"
