// Highway sin on every compiled and supported target: accuracy, bit-identity, comparisons (hwy_test_tmpl.h).
#include <cmath>
extern "C" {
#include "functions/sin/sin_d.h"
}
#define TRIG_NAME "sin"
#define TRIG_HWY_ARRAY hwy_sin_array
#define TRIG_HWY_VARIANT hwy_sin_variant
#define TRIG_HWY_LANES hwy_sin_lanes
#define TRIG_REFL sinl
#define TRIG_STD std::sin
#define TRIG_THRESHOLD 0x1p24
#define TRIG_SVML_ARRAY svml_sin_ha_array
#include "functions/sin/tests/hwy_test_tmpl.h"
