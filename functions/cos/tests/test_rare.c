/* svml_dcos_ha_cout_rare() (C) vs Intel's __svml_dcos_ha_cout_rare_internal, bit for bit. */
#include "../cos_d.h"
#define TRIG_NAME "cos"
#define TRIG_COUT_C svml_dcos_ha_cout_rare
#define TRIG_COUT_ASM __svml_dcos_ha_cout_rare_internal
#include "functions/sin/tests/test_rare_tmpl.h"
