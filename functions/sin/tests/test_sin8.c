/* svml_sin8_ha() (C/AVX-512) vs Intel's __svml_sin8_ha, bit for bit (needs AVX-512, or Intel SDE). */
#include "../sin_d.h"
#define TRIG_NAME "sin"
#define TRIG_C svml_sin8_ha
#define TRIG_ASM __svml_sin8_ha
#define TRIG_THRESHOLD 0x1p24
#include "test_vec_tmpl.h"
