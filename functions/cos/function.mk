# cos, double precision: bit-identical C/AVX-512 translation of svml_z0_cos_d_ha.s and a Google Highway
# version. Shares the large-argument reduction, the constant generator and the test templates with sin
# (functions/sin, functions/sin/trig_common.mk). Targets: cos-all cos-test cos-hwy cos-hashes, and on x86-64 also
# cos-insns cos-avx512-test.

COS     := functions/cos
SIN     := functions/sin
COS_B   := build/cos
COS_ASM := $(SVML_ASM)/svml_z0_cos_d_ha.s
COS_H   := $(COS)/cos_d.h

include functions/sin/trig_common.mk

# --- generated constants (committed; same generator as sin) ---
COS_GEN := $(COS)/svml/svml_cos_d_ha_constants.inc $(TRIG_RED_INC)
$(COS)/svml/svml_cos_d_ha_constants.inc: $(SIN)/tools/gen_svml_constants.py common/tools/extract_tables.py $(COS_ASM)
	python3 -I $< cos $(COS_ASM) > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
COS_OBJS   := $(COS_B)/svml_cos_d_ha.o $(COS_B)/svml_asm.o
COS_TESTS  := $(COS_B)/test_rare $(COS_B)/test_cos8
else
COS_OBJS   := $(COS_B)/svml_cos_d_ha.o
COS_TESTS  :=
endif
COS_PROGS := $(COS_TESTS) $(COS_B)/hwy_test

$(COS_B):
	mkdir -p $@

$(COS_B)/svml_cos_d_ha.o: $(COS)/svml/svml_cos_d_ha.c $(COS_H) $(COS_GEN) $(TRIG_SVML_DEPS) | $(COS_B)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(COS_B)/svml_asm.o: $(COS_ASM) | $(COS_B)
	sed 's/^__svml_dcos_ha_cout_rare_internal:/\t.globl __svml_dcos_ha_cout_rare_internal\n&/' $< > $(COS_B)/svml_asm.s
	$(CC) -c $(COS_B)/svml_asm.s -o $@

# --- tests of the C translation (x86-64) ---
$(COS_B)/test_%: $(COS)/tests/test_%.c $(TRIG_TEST_DEPS) $(COS_H) $(COS_OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(COS_OBJS) -lm -o $@

# --- Highway version ---
$(COS_B)/hwy_cos.o: $(COS)/hwy/hwy_cos.cc $(COS_GEN) $(TRIG_HWY_DEPS) | $(COS_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(COS_B)/extras.o: $(COS)/tests/extras.c $(SIN)/tests/extras_tmpl.h $(SIN)/tests/extras_tmpl_decl.h $(COS_H) | $(COS_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(COS_B)/hwy_test: $(COS)/tests/hwy_test.cc $(TRIG_HWYTEST_DEPS) $(COS_H) $(COS_OBJS) $(COS_B)/hwy_cos.o $(COS_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(COS_B)/hwy_cos.o $(COS_OBJS) $(COS_B)/extras.o $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
cos-all: $(COS_PROGS)
cos-test: $(COS_TESTS) $(COS_B)/hwy_test
	@set -e; for t in $(COS_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(COS_B)/hwy_test 2000000 --no-bench
cos-hwy: $(COS_B)/hwy_test
	./$(COS_B)/hwy_test
cos-hashes: $(COS_B)/hwy_test
	$(RUN) ./$(COS_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += cos-all
TEST_TARGETS   += cos-test
HWY_TARGETS    += cos-hwy
HASHES_TARGETS += cos-hashes

ifeq ($(ARCH),x86_64)
cos-insns: $(COS_B)/svml_cos_d_ha.o $(COS_B)/svml_asm.o
	common/tools/compare_insns.sh $(COS_B)/svml_cos_d_ha.o svml_cos8_ha $(COS_B)/svml_asm.o __svml_cos8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
cos-avx512-test: $(COS_B)/test_cos8 $(COS_B)/hwy_test
	$(RUN) ./$(COS_B)/test_cos8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(COS_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += cos-insns
AVX512_TARGETS += cos-avx512-test
endif

.PHONY: cos-all cos-test cos-hwy cos-hashes cos-insns cos-avx512-test
