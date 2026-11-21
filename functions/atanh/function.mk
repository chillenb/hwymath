# atanh, double precision: bit-identical C/AVX-512 translation of svml_z0_atanh_d_ha.s and a Google Highway
# version that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths
# are relative to the repository root. Targets: atanh-all atanh-test atanh-hwy atanh-hashes atanh-xcheck, and on
# x86-64 also atanh-insns atanh-avx512-test.
# Uses the exact VRCP14PD/VGETEXPPD/VPERMT2PD emulations of common/include/hwymath/hwy/avx512_emul-inl.h (tested
# by common/tests/approx14/test_approx14).

ATANH     := functions/atanh
ATANH_B   := build/atanh
ATANH_ASM := $(SVML_ASM)/svml_z0_atanh_d_ha.s
ATANH_H   := $(ATANH)/atanh_d.h
ATANH_EMUL := $(AVX512_EMUL)

# --- generated constants (committed; regenerated when their generators change) ---
$(ATANH)/svml/svml_atanh_d_ha_constants.inc: $(ATANH)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ATANH_ASM)
	python3 -I $< $(ATANH_ASM) > $@
ATANH_GEN := $(ATANH)/svml/svml_atanh_d_ha_constants.inc

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
ATANH_OBJS   := $(ATANH_B)/svml_atanh_d_ha.o $(ATANH_B)/svml_asm.o
ATANH_EXTRAS := $(ATANH_B)/extras_x86.o
ATANH_TESTS  := $(ATANH_B)/test_rare $(ATANH_B)/test_atanh8
else
ATANH_OBJS   := $(ATANH_B)/svml_atanh_d_ha.o
ATANH_EXTRAS := $(ATANH_B)/extras_arm.o $(ATANH_B)/extras_arm_sve.o
ATANH_TESTS  :=
endif
ATANH_PROGS := $(ATANH_TESTS) $(ATANH_B)/hwy_test

$(ATANH_B):
	mkdir -p $@

# --- C translation ---
$(ATANH_B)/%.o: $(ATANH)/svml/%.c $(ATANH_H) $(ATANH_GEN) | $(ATANH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ATANH_B)/svml_asm.o: $(ATANH_ASM) | $(ATANH_B)
	sed 's/^__svml_datanh_ha_cout_rare_internal:/\t.globl __svml_datanh_ha_cout_rare_internal\n&/' $< > $(ATANH_B)/svml_asm.s
	$(CC) -c $(ATANH_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(ATANH_B)/extras_x86.o: $(ATANH)/tests/extras_x86.c $(ATANH)/tests/extras.h $(ATANH_H) | $(ATANH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ATANH_B)/extras_arm.o: $(ATANH)/tests/extras_arm.c $(ATANH)/tests/extras.h | $(ATANH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
$(ATANH_B)/extras_arm_sve.o: $(ATANH)/tests/extras_arm_sve.c | $(ATANH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
ATANH_TEST_DEPS := $(ATANH)/tests/inputs.h common/include/hwymath/test_util.h $(ATANH_OBJS)
$(ATANH_B)/test_%: $(ATANH)/tests/test_%.c $(ATANH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ATANH_OBJS) -lm -o $@

# --- Highway version ---
# GCC: schedule before register allocation (pressure-aware), so that the two vectors evaluated per loop
# iteration are interleaved; without it one long dependency chain at a time runs (~1.5x slower on AVX2)
HWYMATH_SCHED ?= $(if $(findstring clang,$(shell $(CXX) --version 2>/dev/null | head -1)),,-fschedule-insns -fsched-pressure)
$(ATANH_B)/hwy_atanh.o: $(ATANH)/hwy/hwy_atanh.cc $(ATANH_EMUL) $(ATANH_GEN) | $(ATANH_B)
	$(CXX) $(CXXFLAGS) $(HWYMATH_SCHED) $(INC) -c $< -o $@
$(ATANH_B)/hwy_test: $(ATANH)/tests/hwy_test.cc $(ATANH)/tests/extras.h $(ATANH_TEST_DEPS) $(ATANH_B)/hwy_atanh.o $(ATANH_EXTRAS) common/include/hwymath/accuracy.h
	$(CXX) $(CXXFLAGS) $(INC) $< $(ATANH_B)/hwy_atanh.o $(ATANH_OBJS) $(ATANH_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) -lm -o $@

# --- entry points ---
atanh-all: $(ATANH_PROGS)
atanh-test: $(ATANH_PROGS)
	@set -e; for t in $(ATANH_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ATANH_B)/hwy_test 2000000 --no-bench
atanh-hwy: $(ATANH_B)/hwy_test
	./$(ATANH_B)/hwy_test
atanh-hashes: $(ATANH_B)/hwy_test
	$(RUN) ./$(ATANH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 AVX2, bit for bit (pixi run -e cross make atanh-xcheck)
atanh-xcheck: $(ATANH_GEN)
	XC_SRC="functions/atanh/hwy/hwy_atanh.cc functions/asinh/tests/hwy_xcheck.cc" XC_CXXFLAGS=-DXFN=atanh \
	  XC_VLS="16 32 48 64 256" XC_SVE2_VLS=32 common/tools/xcheck.sh atanh $(PREFIX) $(PREFIX)
ALL_TARGETS    += atanh-all
TEST_TARGETS   += atanh-test
HWY_TARGETS    += atanh-hwy
HASHES_TARGETS += atanh-hashes
XCHECK_TARGETS += atanh-xcheck

ifeq ($(ARCH),x86_64)
atanh-insns: $(ATANH_B)/svml_atanh_d_ha.o $(ATANH_B)/svml_asm.o
	common/tools/compare_insns.sh $(ATANH_B)/svml_atanh_d_ha.o svml_atanh8_ha $(ATANH_B)/svml_asm.o __svml_atanh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly and of every Highway target (AVX-512 ones
# included) vs the translation
atanh-avx512-test: $(ATANH_B)/test_atanh8 $(ATANH_B)/hwy_test
	$(RUN) ./$(ATANH_B)/test_atanh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ATANH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += atanh-insns
AVX512_TARGETS += atanh-avx512-test
endif

.PHONY: atanh-all atanh-test atanh-hwy atanh-hashes atanh-xcheck atanh-insns atanh-avx512-test
