# Make variables and rules shared by sin, cos and tan (functions/{sin,cos,tan}/function.mk include this
# file; the top-level Makefile includes the function.mk files in alphabetical order, so cos comes first).
ifndef TRIG_COMMON_MK
TRIG_COMMON_MK := 1

TRIG := functions/sin
# large-argument reduction table and constants, generated from the sin assembly (identical in all three)
TRIG_RED_INC := $(TRIG)/svml/svml_trig_reduction_data.inc
$(TRIG_RED_INC): $(TRIG)/tools/gen_svml_constants.py common/tools/extract_tables.py $(SVML_ASM)/svml_z0_sin_d_ha.s
	python3 -I $< reduction $(SVML_ASM)/svml_z0_sin_d_ha.s > $@

TRIG_SVML_DEPS    := $(TRIG_RED_INC) $(TRIG)/svml/svml_trig_reduction.h
TRIG_HWY_DEPS     := $(TRIG_RED_INC) $(TRIG)/hwy/trig_hwy-inl.h
TRIG_TEST_DEPS    := $(TRIG)/tests/inputs.h $(TRIG)/tests/test_rare_tmpl.h $(TRIG)/tests/test_vec_tmpl.h \
                     common/include/hwymath/test_util.h
TRIG_HWYTEST_DEPS := $(TRIG)/tests/inputs.h $(TRIG)/tests/hwy_test_tmpl.h common/include/hwymath/test_util.h \
                     common/include/hwymath/accuracy.h
endif
