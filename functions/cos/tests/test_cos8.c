/* svml_cos8_ha() (C/AVX-512) vs Intel's __svml_cos8_ha, bit for bit (needs AVX-512, or Intel SDE). */
#include "../cos_d.h"
#define TRIG_NAME "cos"
#define TRIG_C svml_cos8_ha
#define TRIG_ASM __svml_cos8_ha
#define TRIG_THRESHOLD 0x1p24
#include "functions/sin/tests/test_vec_tmpl.h"
