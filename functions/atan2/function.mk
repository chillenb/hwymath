# atan2, double precision: bit-identical C/AVX-512 translation of svml_z0_atan2_d_ha.s and a Google Highway
# version. Included by the top-level Makefile; paths are relative to the repository root.
# Targets: atan2-all atan2-test atan2-hwy atan2-hashes atan2-xcheck, and on x86-64 also atan2-insns atan2-avx512-test.
# The scalar callout is translated by functions/asin/tools/asm2c.py (machine model: functions/asin/svml/asm_model.h).

ATAN2     := functions/atan2
ATAN2_B   := build/atan2
ATAN2_ASM := $(SVML_ASM)/svml_z0_atan2_d_ha.s
ATAN2_H   := $(ATAN2)/atan2_d.h
ATAN2_ASM2C   := functions/asin/tools/asm2c.py
ATAN2_MODEL_H := functions/asin/svml/asm_model.h

# --- generated sources (committed; regenerated when their generators change) ---
ATAN2_GEN := $(ATAN2)/svml/svml_atan2_d_ha_constants.inc $(ATAN2)/svml/svml_datan2_ha_cout_rare.inc
$(ATAN2)/svml/svml_atan2_d_ha_constants.inc: $(ATAN2)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ATAN2_ASM)
	python3 -I $< $(ATAN2_ASM) > $@
$(ATAN2)/svml/svml_datan2_ha_cout_rare.inc: $(ATAN2_ASM2C) common/tools/extract_tables.py $(ATAN2_ASM)
	python3 -I $(ATAN2_ASM2C) $(ATAN2_ASM) __svml_datan2_ha_cout_rare_internal svml_datan2_ha_cout_rare 3 > $@

# --- what gets built on this architecture ---
ATAN2_OBJS := $(ATAN2_B)/svml_atan2_d_ha.o
ifeq ($(ARCH),x86_64)
ATAN2_OBJS  += $(ATAN2_B)/svml_asm.o
ATAN2_TESTS := $(ATAN2_B)/test_rare $(ATAN2_B)/test_atan28
else
ATAN2_TESTS :=
endif
ATAN2_PROGS := $(ATAN2_TESTS) $(ATAN2_B)/hwy_test

$(ATAN2_B):
	mkdir -p $@

$(ATAN2_B)/svml_atan2_d_ha.o: $(ATAN2)/svml/svml_atan2_d_ha.c $(ATAN2_H) $(ATAN2_MODEL_H) $(ATAN2_GEN) | $(ATAN2_B)
	$(CC) $(CFLAGS) -I$(dir $(ATAN2_MODEL_H)) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ATAN2_B)/svml_asm.o: $(ATAN2_ASM) | $(ATAN2_B)
	sed 's/^__svml_datan2_ha_cout_rare_internal:/\t.globl __svml_datan2_ha_cout_rare_internal\n&/' $< > $(ATAN2_B)/svml_asm.s
	$(CC) -c $(ATAN2_B)/svml_asm.s -o $@

ATAN2_TEST_DEPS := $(ATAN2)/tests/inputs.h common/include/hwymath/test_util.h $(ATAN2_OBJS)
$(ATAN2_B)/test_%: $(ATAN2)/tests/test_%.c $(ATAN2_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ATAN2_OBJS) -lm -o $@

# --- Highway version ---
$(ATAN2_B)/hwy_atan2.o: $(ATAN2)/hwy/hwy_atan2.cc $(ATAN2)/hwy/hwy_atan2_constants.h $(ATAN2_GEN) $(AVX512_EMUL) | $(ATAN2_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(ATAN2_B)/extras.o: $(ATAN2)/tests/extras.c $(ATAN2)/tests/extras.h | $(ATAN2_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ATAN2_B)/hwy_test: $(ATAN2)/tests/hwy_test.cc $(ATAN2)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(ATAN2_TEST_DEPS) $(ATAN2_B)/hwy_atan2.o $(ATAN2_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(ATAN2_B)/hwy_atan2.o $(ATAN2_OBJS) $(ATAN2_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

# --- entry points ---
atan2-all: $(ATAN2_PROGS)
atan2-test: $(ATAN2_TESTS) $(ATAN2_B)/hwy_test
	@set -e; for t in $(ATAN2_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ATAN2_B)/hwy_test 2000000 --no-bench
atan2-hwy: $(ATAN2_B)/hwy_test
	./$(ATAN2_B)/hwy_test
atan2-hashes: $(ATAN2_B)/hwy_test
	$(RUN) ./$(ATAN2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += atan2-all
TEST_TARGETS   += atan2-test
HWY_TARGETS    += atan2-hwy
HASHES_TARGETS += atan2-hashes
# AArch64 NEON/SVE/SVE2 (QEMU, SVE at 128-2048 bits) vs x86-64, bit for bit: pixi run -e cross make atan2-xcheck
atan2-xcheck: $(ATAN2_GEN)
	XC_SRC="functions/atan2/hwy/hwy_atan2.cc functions/asin/tools/xcheck_invtrig.cc" XC_CFLAGS=-Ifunctions/asin/svml \
	  XC_CXXFLAGS="-Ifunctions/atan2/tests -DXC_F=atan2 -DXC_ARGS=2" XC_VLS="16 32 64 128 256" \
	  common/tools/xcheck.sh atan2 $(PREFIX) $(PREFIX) 200000
XCHECK_TARGETS += atan2-xcheck

ifeq ($(ARCH),x86_64)
atan2-insns: $(ATAN2_B)/svml_atan2_d_ha.o $(ATAN2_B)/svml_asm.o
	common/tools/compare_insns.sh $(ATAN2_B)/svml_atan2_d_ha.o svml_atan28_ha $(ATAN2_B)/svml_asm.o __svml_atan28_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
atan2-avx512-test: $(ATAN2_B)/test_atan28 $(ATAN2_B)/hwy_test
	$(RUN) ./$(ATAN2_B)/test_atan28 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ATAN2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += atan2-insns
AVX512_TARGETS += atan2-avx512-test
endif

.PHONY: atan2-all atan2-test atan2-hwy atan2-hashes atan2-xcheck atan2-insns atan2-avx512-test
