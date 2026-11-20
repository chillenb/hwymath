/* svml_tan8_ha() (C/AVX-512) vs Intel's __svml_tan8_ha, bit for bit (needs AVX-512, or Intel SDE). */
#include "../tan_d.h"
#define TRIG_NAME "tan"
#define TRIG_C svml_tan8_ha
#define TRIG_ASM __svml_tan8_ha
#define TRIG_THRESHOLD 0x1.000000e4db24cp+16
#include "functions/sin/tests/test_vec_tmpl.h"
