# tan, double precision: bit-identical C/AVX-512 translation of svml_z0_tan_d_ha.s and a Google Highway
# version (with an exact emulation of vrcp14pd). Shares the large-argument reduction and the test
# templates with sin (functions/sin, functions/sin/trig_common.mk). Targets: tan-all tan-test tan-hwy
# tan-hashes tan-xcheck, and on x86-64 also tan-insns tan-avx512-test.

TAN     := functions/tan
SIN     := functions/sin
TAN_B   := build/tan
TAN_ASM := $(SVML_ASM)/svml_z0_tan_d_ha.s
TAN_H   := $(TAN)/tan_d.h

include functions/sin/trig_common.mk

# --- generated constants (committed; regenerated when their generator changes) ---
TAN_GEN := $(TAN)/svml/svml_tan_d_ha_constants.inc $(TRIG_RED_INC)
$(TAN)/svml/svml_tan_d_ha_constants.inc: $(TAN)/tools/gen_svml_constants.py common/tools/extract_tables.py $(TAN_ASM)
	python3 -I $< $(TAN_ASM) > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
TAN_OBJS   := $(TAN_B)/svml_tan_d_ha.o $(TAN_B)/svml_asm.o
TAN_TESTS  := $(TAN_B)/test_rare $(TAN_B)/test_rcp14 $(TAN_B)/test_tan8
else
TAN_OBJS   := $(TAN_B)/svml_tan_d_ha.o
TAN_TESTS  :=
endif
TAN_PROGS := $(TAN_TESTS) $(TAN_B)/hwy_test

$(TAN_B):
	mkdir -p $@

$(TAN_B)/svml_tan_d_ha.o: $(TAN)/svml/svml_tan_d_ha.c $(TAN_H) $(TAN_GEN) $(TRIG_SVML_DEPS) | $(TAN_B)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly, and a global
# label __svml_tan8_ha_eval at the start of the evaluation part (.LBL_1_2) for make tan-insns
$(TAN_B)/svml_asm.o: $(TAN_ASM) | $(TAN_B)
	sed -e 's/^__svml_dtan_ha_cout_rare_internal:/\t.globl __svml_dtan_ha_cout_rare_internal\n&/' \
	    -e 's/^\.LBL_1_2:/&\n\t.globl __svml_tan8_ha_eval\n__svml_tan8_ha_eval:/' $< > $(TAN_B)/svml_asm.s
	$(CC) -c $(TAN_B)/svml_asm.s -o $@

# --- tests of the C translation (x86-64) ---
$(TAN_B)/test_%: $(TAN)/tests/test_%.c $(TRIG_TEST_DEPS) common/include/hwymath/rcp14_ref.h $(TAN)/tests/rcp14_inputs.h $(TAN_H) $(TAN_OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(TAN_OBJS) -lm -o $@

# --- Highway version ---
$(TAN_B)/hwy_tan.o: $(TAN)/hwy/hwy_tan.cc $(TAN)/hwy/rcp14_hwy-inl.h common/include/hwymath/rcp14_ref.h $(TAN_GEN) $(TRIG_HWY_DEPS) | $(TAN_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(TAN_B)/extras.o: $(TAN)/tests/extras.c $(SIN)/tests/extras_tmpl.h $(SIN)/tests/extras_tmpl_decl.h $(TAN_H) | $(TAN_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(TAN_B)/hwy_test: $(TAN)/tests/hwy_test.cc $(TRIG_HWYTEST_DEPS) common/include/hwymath/rcp14_ref.h $(TAN)/tests/rcp14_inputs.h $(TAN_H) $(TAN_OBJS) $(TAN_B)/hwy_tan.o $(TAN_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(TAN_B)/hwy_tan.o $(TAN_OBJS) $(TAN_B)/extras.o $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
tan-all: $(TAN_PROGS)
tan-test: $(TAN_TESTS) $(TAN_B)/hwy_test
	@set -e; for t in $(TAN_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(TAN_B)/hwy_test 2000000 --no-bench
tan-hwy: $(TAN_B)/hwy_test
	./$(TAN_B)/hwy_test
tan-hashes: $(TAN_B)/hwy_test
	$(RUN) ./$(TAN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += tan-all
TEST_TARGETS   += tan-test
HWY_TARGETS    += tan-hwy
HASHES_TARGETS += tan-hashes
# Highway tan on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (run in the cross environment:
# pixi run -e cross make tan-xcheck)
tan-xcheck: $(TAN_GEN)
	XC_TARGETS="neon:-march=armv8-a neon_aes:-march=armv8-a+crypto sve:-march=armv8.2-a+sve sve2:-march=armv9-a+sve2" \
	  common/tools/xcheck.sh tan $(PREFIX) $(PREFIX)
XCHECK_TARGETS += tan-xcheck

ifeq ($(ARCH),x86_64)
# two parts: the reduction up to the branch to the large path, and the evaluation (.LBL_1_2)
tan-insns: $(TAN_B)/svml_tan_d_ha.o $(TAN_B)/svml_asm.o
	common/tools/compare_insns.sh $(TAN_B)/svml_tan_d_ha.o svml_tan8_ha $(TAN_B)/svml_asm.o __svml_tan8_ha
	common/tools/compare_insns.sh $(TAN_B)/svml_tan_d_ha.o svml_tan8_ha_eval $(TAN_B)/svml_asm.o __svml_tan8_ha_eval
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
tan-avx512-test: $(TAN_B)/test_rcp14 $(TAN_B)/test_tan8 $(TAN_B)/hwy_test
	$(RUN) ./$(TAN_B)/test_rcp14 || [ $$? -eq 77 ]
	$(RUN) ./$(TAN_B)/test_tan8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(TAN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += tan-insns
AVX512_TARGETS += tan-avx512-test
endif

.PHONY: tan-all tan-test tan-hwy tan-hashes tan-xcheck tan-insns tan-avx512-test
