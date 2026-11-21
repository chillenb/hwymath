# asinh, double precision: bit-identical C/AVX-512 translation of svml_z0_asinh_d_ha.s and a Google Highway
# version that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths
# are relative to the repository root. Targets: asinh-all asinh-test asinh-hwy asinh-hashes asinh-xcheck, and on
# x86-64 also asinh-insns asinh-avx512-test.
# Uses the exact AVX-512 emulations of common/include/hwymath/hwy/avx512_emul-inl.h (tested by
# common/tests/approx14/test_approx14).

ASINH     := functions/asinh
ASINH_B   := build/asinh
ASINH_ASM := $(SVML_ASM)/svml_z0_asinh_d_ha.s
ASINH_H   := $(ASINH)/asinh_d.h

# --- generated constants (committed; regenerated when their generators change) ---
$(ASINH)/svml/svml_asinh_d_ha_constants.inc: $(ASINH)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ASINH_ASM)
	python3 -I $< $(ASINH_ASM) asinh > $@
ASINH_GEN := $(ASINH)/svml/svml_asinh_d_ha_constants.inc

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
ASINH_OBJS   := $(ASINH_B)/svml_asinh_d_ha.o $(ASINH_B)/svml_asm.o
ASINH_EXTRAS := $(ASINH_B)/extras_x86.o
ASINH_TESTS  := $(ASINH_B)/test_rare $(ASINH_B)/test_asinh8
else
ASINH_OBJS   := $(ASINH_B)/svml_asinh_d_ha.o
ASINH_EXTRAS := $(ASINH_B)/extras_arm.o $(ASINH_B)/extras_arm_sve.o
ASINH_TESTS  :=
endif
ASINH_PROGS := $(ASINH_TESTS) $(ASINH_B)/hwy_test

$(ASINH_B):
	mkdir -p $@

# --- C translation ---
$(ASINH_B)/%.o: $(ASINH)/svml/%.c $(ASINH_H) $(ASINH_GEN) | $(ASINH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ASINH_B)/svml_asm.o: $(ASINH_ASM) | $(ASINH_B)
	sed 's/^__svml_dasinh_ha_cout_rare_internal:/\t.globl __svml_dasinh_ha_cout_rare_internal\n&/' $< > $(ASINH_B)/svml_asm.s
	$(CC) -c $(ASINH_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(ASINH_B)/extras_x86.o: $(ASINH)/tests/extras_x86.c $(ASINH)/tests/extras.h $(ASINH_H) | $(ASINH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ASINH_B)/extras_arm.o: $(ASINH)/tests/extras_arm.c $(ASINH)/tests/extras.h | $(ASINH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
$(ASINH_B)/extras_arm_sve.o: $(ASINH)/tests/extras_arm_sve.c | $(ASINH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
ASINH_TEST_DEPS := $(ASINH)/tests/inputs.h common/include/hwymath/test_util.h $(ASINH_OBJS)
$(ASINH_B)/test_%: $(ASINH)/tests/test_%.c $(ASINH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ASINH_OBJS) -lm -o $@

# --- Highway version ---
# GCC: schedule before register allocation (pressure-aware), so that the two vectors evaluated per loop
# iteration are interleaved; without it one long dependency chain at a time runs (~1.5x slower on AVX2)
HWYMATH_SCHED ?= $(if $(findstring clang,$(shell $(CXX) --version 2>/dev/null | head -1)),,-fschedule-insns -fsched-pressure)
$(ASINH_B)/hwy_asinh.o: $(ASINH)/hwy/hwy_asinh.cc $(AVX512_EMUL) $(ASINH_GEN) | $(ASINH_B)
	$(CXX) $(CXXFLAGS) $(HWYMATH_SCHED) $(INC) -c $< -o $@
$(ASINH_B)/hwy_test: $(ASINH)/tests/hwy_test.cc $(ASINH)/tests/extras.h $(ASINH_TEST_DEPS) $(ASINH_B)/hwy_asinh.o $(ASINH_EXTRAS) common/include/hwymath/accuracy.h
	$(CXX) $(CXXFLAGS) $(INC) $< $(ASINH_B)/hwy_asinh.o $(ASINH_OBJS) $(ASINH_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) -lm -o $@

# --- entry points ---
asinh-all: $(ASINH_PROGS)
asinh-test: $(ASINH_PROGS)
	@set -e; for t in $(ASINH_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ASINH_B)/hwy_test 2000000 --no-bench
asinh-hwy: $(ASINH_B)/hwy_test
	./$(ASINH_B)/hwy_test
asinh-hashes: $(ASINH_B)/hwy_test
	$(RUN) ./$(ASINH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 AVX2, bit for bit (pixi run -e cross make asinh-xcheck)
asinh-xcheck: $(ASINH_GEN)
	XC_SRC="functions/asinh/hwy/hwy_asinh.cc functions/asinh/tests/hwy_xcheck.cc" XC_CXXFLAGS=-DXFN=asinh \
	  XC_VLS="16 32 48 64 256" XC_SVE2_VLS=32 common/tools/xcheck.sh asinh $(PREFIX) $(PREFIX)
ALL_TARGETS    += asinh-all
TEST_TARGETS   += asinh-test
HWY_TARGETS    += asinh-hwy
HASHES_TARGETS += asinh-hashes
XCHECK_TARGETS += asinh-xcheck

ifeq ($(ARCH),x86_64)
asinh-insns: $(ASINH_B)/svml_asinh_d_ha.o $(ASINH_B)/svml_asm.o
	common/tools/compare_insns.sh $(ASINH_B)/svml_asinh_d_ha.o svml_asinh8_ha $(ASINH_B)/svml_asm.o __svml_asinh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly and of every Highway target (AVX-512 ones included) vs the translation
asinh-avx512-test: $(ASINH_B)/test_asinh8 $(ASINH_B)/hwy_test
	$(RUN) ./$(ASINH_B)/test_asinh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ASINH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += asinh-insns
AVX512_TARGETS += asinh-avx512-test
endif

.PHONY: asinh-all asinh-test asinh-hwy asinh-hashes asinh-xcheck asinh-insns asinh-avx512-test
