// Highway cos on every compiled and supported target: accuracy, bit-identity, comparisons (hwy_test_tmpl.h).
#include <cmath>
extern "C" {
#include "functions/cos/cos_d.h"
}
#define TRIG_NAME "cos"
#define TRIG_HWY_ARRAY hwy_cos_array
#define TRIG_HWY_VARIANT hwy_cos_variant
#define TRIG_HWY_LANES hwy_cos_lanes
#define TRIG_REFL cosl
#define TRIG_STD std::cos
#define TRIG_THRESHOLD 0x1p24
#define TRIG_SVML_ARRAY svml_cos_ha_array
#include "functions/sin/tests/hwy_test_tmpl.h"
