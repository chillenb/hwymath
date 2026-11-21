# cbrt, double precision: bit-identical C/AVX-512 translation of svml_z0_cbrt_d_ha.s and a Google
# Highway version of the same algorithm. Included by the top-level Makefile; paths are relative to the
# repository root. Targets: cbrt-all cbrt-test cbrt-hwy cbrt-hashes, and on x86-64 also cbrt-insns
# cbrt-avx512-test.

CBRT     := functions/cbrt
CBRT_B   := build/cbrt
CBRT_ASM := $(SVML_ASM)/svml_z0_cbrt_d_ha.s
CBRT_H   := $(CBRT)/cbrt_d.h

# --- generated constants (committed; regenerated when the generator changes) ---
CBRT_GEN := $(CBRT)/svml/svml_cbrt_d_ha_constants.inc
$(CBRT)/svml/svml_cbrt_d_ha_constants.inc: $(CBRT)/tools/gen_svml_constants.py common/tools/extract_tables.py $(CBRT_ASM)
	python3 -I $< $(CBRT_ASM) > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
CBRT_OBJS   := $(CBRT_B)/svml_cbrt_d_ha.o $(CBRT_B)/svml_asm.o
CBRT_EXTRAS := $(CBRT_B)/extras_x86.o
CBRT_TESTS  := $(CBRT_B)/test_rare $(CBRT_B)/test_cbrt8 $(CBRT_B)/test_rcp14
else
CBRT_OBJS   := $(CBRT_B)/svml_cbrt_d_ha.o
CBRT_EXTRAS := $(CBRT_B)/extras_arm.o $(CBRT_B)/extras_arm_sve.o
CBRT_TESTS  :=
endif
CBRT_PROGS := $(CBRT_TESTS) $(CBRT_B)/hwy_test

$(CBRT_B):
	mkdir -p $@

# --- C translation ---
$(CBRT_B)/svml_cbrt_d_ha.o: $(CBRT)/svml/svml_cbrt_d_ha.c $(CBRT_H) $(CBRT_GEN) | $(CBRT_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(CBRT_B)/svml_asm.o: $(CBRT_ASM) | $(CBRT_B)
	sed 's/^__svml_dcbrt_ha_cout_rare_internal:/\t.globl __svml_dcbrt_ha_cout_rare_internal\n&/' $< > $(CBRT_B)/svml_asm.s
	$(CC) -c $(CBRT_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(CBRT_B)/extras_x86.o: $(CBRT)/tests/extras_x86.c $(CBRT)/tests/extras.h $(CBRT_H) | $(CBRT_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(CBRT_B)/extras_arm.o: $(CBRT)/tests/extras_arm.c $(CBRT)/tests/extras.h $(CBRT_H) | $(CBRT_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
$(CBRT_B)/extras_arm_sve.o: $(CBRT)/tests/extras_arm_sve.c | $(CBRT_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
CBRT_TEST_DEPS := $(CBRT)/tests/inputs.h $(CBRT)/rcp14_model.h common/include/hwymath/test_util.h $(CBRT_OBJS)
$(CBRT_B)/test_rcp14: $(CBRT)/tests/test_rcp14.c $(CBRT)/rcp14_model.h common/include/hwymath/test_util.h | $(CBRT_B)
	$(CC) $(CFLAGS) $(INC) $< -o $@
$(CBRT_B)/test_%: $(CBRT)/tests/test_%.c $(CBRT_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(CBRT_OBJS) -lm -o $@

# --- Highway version ---
$(CBRT_B)/hwy_cbrt.o: $(CBRT)/hwy/hwy_cbrt.cc common/include/hwymath/hwy/avx512_native-inl.h $(CBRT)/hwy/hwy_cbrt_constants.h $(CBRT)/rcp14_model.h $(CBRT_GEN) | $(CBRT_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(CBRT_B)/hwy_test: $(CBRT)/tests/hwy_test.cc $(CBRT)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(CBRT_TEST_DEPS) $(CBRT_B)/hwy_cbrt.o $(CBRT_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(CBRT_B)/hwy_cbrt.o $(CBRT_OBJS) $(CBRT_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(if $(filter x86_64,$(ARCH)),$(INTEL_LIB)) -lm -o $@

# --- entry points ---
cbrt-all: $(CBRT_PROGS)
cbrt-test: $(CBRT_TESTS) $(CBRT_B)/hwy_test
	@set -e; for t in $(CBRT_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(CBRT_B)/hwy_test 2000000 --no-bench
cbrt-hwy: $(CBRT_B)/hwy_test
	./$(CBRT_B)/hwy_test
cbrt-hashes: $(CBRT_B)/hwy_test
	$(RUN) ./$(CBRT_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += cbrt-all
TEST_TARGETS   += cbrt-test
HWY_TARGETS    += cbrt-hwy
HASHES_TARGETS += cbrt-hashes

ifeq ($(ARCH),x86_64)
cbrt-insns: $(CBRT_B)/svml_cbrt_d_ha.o $(CBRT_B)/svml_asm.o
	common/tools/compare_insns.sh $(CBRT_B)/svml_cbrt_d_ha.o svml_cbrt8_ha $(CBRT_B)/svml_asm.o __svml_cbrt8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, the rcp14 model, Highway AVX-512 targets
# (which hwy_test also compares bit for bit with the translation)
cbrt-avx512-test: $(CBRT_B)/test_cbrt8 $(CBRT_B)/test_rcp14 $(CBRT_B)/hwy_test
	$(RUN) ./$(CBRT_B)/test_cbrt8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(CBRT_B)/test_rcp14 || [ $$? -eq 77 ]
	$(RUN) ./$(CBRT_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += cbrt-insns
AVX512_TARGETS += cbrt-avx512-test
endif

.PHONY: cbrt-all cbrt-test cbrt-hwy cbrt-hashes cbrt-insns cbrt-avx512-test
