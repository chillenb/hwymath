/* svml_dtan_ha_cout_rare() (C) vs Intel's __svml_dtan_ha_cout_rare_internal, bit for bit. */
#include "../tan_d.h"
#define TRIG_NAME "tan"
#define TRIG_COUT_C svml_dtan_ha_cout_rare
#define TRIG_COUT_ASM __svml_dtan_ha_cout_rare_internal
#include "functions/sin/tests/test_rare_tmpl.h"
