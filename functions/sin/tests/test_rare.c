/* svml_dsin_ha_cout_rare() (C) vs Intel's __svml_dsin_ha_cout_rare_internal, bit for bit. */
#include "../sin_d.h"
#define TRIG_NAME "sin"
#define TRIG_COUT_C svml_dsin_ha_cout_rare
#define TRIG_COUT_ASM __svml_dsin_ha_cout_rare_internal
#include "test_rare_tmpl.h"
