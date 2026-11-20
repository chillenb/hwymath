# asin, double precision: bit-identical C/AVX-512 translation of svml_z0_asin_d_ha.s and a Google Highway
# version. Included by the top-level Makefile; paths are relative to the repository root.
# Targets: asin-all asin-test asin-hwy asin-hashes asin-xcheck, and on x86-64 also asin-insns asin-avx512-test.
# functions/asin/tools/asm2c.py and functions/asin/svml/asm_model.h are also used by acos, atan, atan2.

ASIN     := functions/asin
ASIN_B   := build/asin
ASIN_ASM := $(SVML_ASM)/svml_z0_asin_d_ha.s
ASIN_H   := $(ASIN)/asin_d.h
ASIN_ASM2C  := functions/asin/tools/asm2c.py
ASIN_MODEL_H := functions/asin/svml/asm_model.h

# --- generated sources (committed; regenerated when their generators change) ---
ASIN_GEN := $(ASIN)/svml/svml_asin_d_ha_constants.inc $(ASIN)/svml/svml_dasin_ha_cout_rare.inc
$(ASIN)/svml/svml_asin_d_ha_constants.inc: $(ASIN)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ASIN_ASM)
	python3 -I $< $(ASIN_ASM) > $@
$(ASIN)/svml/svml_dasin_ha_cout_rare.inc: $(ASIN_ASM2C) common/tools/extract_tables.py $(ASIN_ASM)
	python3 -I $(ASIN_ASM2C) $(ASIN_ASM) __svml_dasin_ha_cout_rare_internal svml_dasin_ha_cout_rare 2 > $@

# --- what gets built on this architecture ---
ASIN_OBJS := $(ASIN_B)/svml_asin_d_ha.o
ifeq ($(ARCH),x86_64)
ASIN_OBJS  += $(ASIN_B)/svml_asm.o
ASIN_TESTS := $(ASIN_B)/test_rare $(ASIN_B)/test_asin8
else
ASIN_TESTS :=
endif
ASIN_PROGS := $(ASIN_TESTS) $(ASIN_B)/test_approx14 $(ASIN_B)/hwy_test

$(ASIN_B):
	mkdir -p $@

$(ASIN_B)/svml_asin_d_ha.o: $(ASIN)/svml/svml_asin_d_ha.c $(ASIN_H) $(ASIN_MODEL_H) $(ASIN_GEN) | $(ASIN_B)
	$(CC) $(CFLAGS) -I$(dir $(ASIN_MODEL_H)) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ASIN_B)/svml_asm.o: $(ASIN_ASM) | $(ASIN_B)
	sed 's/^__svml_dasin_ha_cout_rare_internal:/\t.globl __svml_dasin_ha_cout_rare_internal\n&/' $< > $(ASIN_B)/svml_asm.s
	$(CC) -c $(ASIN_B)/svml_asm.s -o $@

ASIN_TEST_DEPS := $(ASIN)/tests/inputs.h common/include/hwymath/test_util.h $(ASIN_OBJS)
$(ASIN_B)/test_%: $(ASIN)/tests/test_%.c $(ASIN_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ASIN_OBJS) -lm -o $@

# --- exact emulation of vrcp14pd/vrsqrt14pd (used by asin, acos, atan, atan2) and its exhaustive test ---
functions/asin/hwy/approx14_tables.h: functions/asin/tools/gen_approx14.py
	python3 -I $< > $@
$(ASIN_B)/hwy_approx14.o: $(ASIN)/tests/hwy_approx14.cc functions/asin/hwy/approx14-inl.h functions/asin/hwy/approx14_tables.h | $(ASIN_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(ASIN_B)/test_approx14: $(ASIN)/tests/test_approx14.cc $(ASIN)/tests/approx14_ref.h $(ASIN_B)/hwy_approx14.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(ASIN_B)/hwy_approx14.o $(LIBPATH) -lhwy -lm -o $@

# --- Highway version ---
$(ASIN_B)/hwy_asin.o: $(ASIN)/hwy/hwy_asin.cc $(ASIN)/hwy/hwy_asin_constants.h $(ASIN_GEN) functions/asin/hwy/approx14-inl.h functions/asin/hwy/approx14_tables.h | $(ASIN_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(ASIN_B)/extras.o: $(ASIN)/tests/extras.c $(ASIN)/tests/extras.h | $(ASIN_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ASIN_B)/hwy_test: $(ASIN)/tests/hwy_test.cc $(ASIN)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(ASIN_TEST_DEPS) $(ASIN_B)/hwy_asin.o $(ASIN_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(ASIN_B)/hwy_asin.o $(ASIN_OBJS) $(ASIN_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

# --- entry points ---
asin-all: $(ASIN_PROGS)
asin-test: $(ASIN_TESTS) build/asin/test_approx14 $(ASIN_B)/hwy_test
	@set -e; for t in $(ASIN_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./build/asin/test_approx14
	$(RUN) ./$(ASIN_B)/hwy_test 2000000 --no-bench
asin-hwy: $(ASIN_B)/hwy_test
	./$(ASIN_B)/hwy_test
asin-hashes: $(ASIN_B)/hwy_test
	$(RUN) ./$(ASIN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += asin-all
TEST_TARGETS   += asin-test
HWY_TARGETS    += asin-hwy
HASHES_TARGETS += asin-hashes
# AArch64 NEON/SVE/SVE2 (QEMU, SVE at 128-2048 bits) vs x86-64, bit for bit: pixi run -e cross make asin-xcheck
asin-xcheck: $(ASIN_GEN)
	functions/asin/tools/xcheck_invtrig.sh asin $(PREFIX) $(PREFIX)
XCHECK_TARGETS += asin-xcheck

ifeq ($(ARCH),x86_64)
asin-insns: $(ASIN_B)/svml_asin_d_ha.o $(ASIN_B)/svml_asm.o
	common/tools/compare_insns.sh $(ASIN_B)/svml_asin_d_ha.o svml_asin8_ha $(ASIN_B)/svml_asm.o __svml_asin8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
asin-avx512-test: $(ASIN_B)/test_asin8 $(ASIN_B)/hwy_test build/asin/test_approx14
	$(RUN) ./build/asin/test_approx14
	$(RUN) ./$(ASIN_B)/test_asin8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ASIN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += asin-insns
AVX512_TARGETS += asin-avx512-test
endif

.PHONY: asin-all asin-test asin-hwy asin-hashes asin-xcheck asin-insns asin-avx512-test
