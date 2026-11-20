# acosh, double precision: bit-identical C/AVX-512 translation of svml_z0_acosh_d_ha.s and a Google Highway
# version that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths
# are relative to the repository root. Targets: acosh-all acosh-test acosh-hwy acosh-hashes acosh-xcheck, and on
# x86-64 also acosh-insns acosh-avx512-test.
# Uses the exact VRCP14PD/VRSQRT14PD/VGETEXPPD/VPERMT2PD emulations of functions/asinh/hwy/avx512_emul-inl.h
# (tested by asinh's test_approx14) and asinh's constant generator (same data layout).

ACOSH     := functions/acosh
ACOSH_B   := build/acosh
ACOSH_ASM := $(SVML_ASM)/svml_z0_acosh_d_ha.s
ACOSH_H   := $(ACOSH)/acosh_d.h
ACOSH_EMUL := functions/asinh/hwy/avx512_emul-inl.h functions/asinh/hwy/approx14_tables.h

# --- generated constants (committed; regenerated when their generators change) ---
$(ACOSH)/svml/svml_acosh_d_ha_constants.inc: $(ACOSH)/tools/gen_svml_constants.py functions/asinh/tools/gen_svml_constants.py common/tools/extract_tables.py $(ACOSH_ASM)
	python3 -I $< $(ACOSH_ASM) > $@
ACOSH_GEN := $(ACOSH)/svml/svml_acosh_d_ha_constants.inc

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
ACOSH_OBJS   := $(ACOSH_B)/svml_acosh_d_ha.o $(ACOSH_B)/svml_asm.o
ACOSH_EXTRAS := $(ACOSH_B)/extras_x86.o
ACOSH_TESTS  := $(ACOSH_B)/test_rare $(ACOSH_B)/test_acosh8
else
ACOSH_OBJS   := $(ACOSH_B)/svml_acosh_d_ha.o
ACOSH_EXTRAS := $(ACOSH_B)/extras_arm.o $(ACOSH_B)/extras_arm_sve.o
ACOSH_TESTS  :=
endif
ACOSH_PROGS := $(ACOSH_TESTS) $(ACOSH_B)/hwy_test

$(ACOSH_B):
	mkdir -p $@

# --- C translation ---
$(ACOSH_B)/%.o: $(ACOSH)/svml/%.c $(ACOSH_H) $(ACOSH_GEN) | $(ACOSH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ACOSH_B)/svml_asm.o: $(ACOSH_ASM) | $(ACOSH_B)
	sed 's/^__svml_dacosh_ha_cout_rare_internal:/\t.globl __svml_dacosh_ha_cout_rare_internal\n&/' $< > $(ACOSH_B)/svml_asm.s
	$(CC) -c $(ACOSH_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(ACOSH_B)/extras_x86.o: $(ACOSH)/tests/extras_x86.c $(ACOSH)/tests/extras.h $(ACOSH_H) | $(ACOSH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ACOSH_B)/extras_arm.o: $(ACOSH)/tests/extras_arm.c $(ACOSH)/tests/extras.h | $(ACOSH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
$(ACOSH_B)/extras_arm_sve.o: $(ACOSH)/tests/extras_arm_sve.c | $(ACOSH_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
ACOSH_TEST_DEPS := $(ACOSH)/tests/inputs.h common/include/hwymath/test_util.h $(ACOSH_OBJS)
$(ACOSH_B)/test_%: $(ACOSH)/tests/test_%.c $(ACOSH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ACOSH_OBJS) -lm -o $@

# --- Highway version ---
# GCC: schedule before register allocation (pressure-aware), so that the two vectors evaluated per loop
# iteration are interleaved; without it one long dependency chain at a time runs (~1.5x slower on AVX2)
HWYMATH_SCHED ?= $(if $(findstring clang,$(shell $(CXX) --version 2>/dev/null | head -1)),,-fschedule-insns -fsched-pressure)
$(ACOSH_B)/hwy_acosh.o: $(ACOSH)/hwy/hwy_acosh.cc $(ACOSH_EMUL) $(ACOSH_GEN) | $(ACOSH_B)
	$(CXX) $(CXXFLAGS) $(HWYMATH_SCHED) $(INC) -c $< -o $@
$(ACOSH_B)/hwy_test: $(ACOSH)/tests/hwy_test.cc $(ACOSH)/tests/extras.h $(ACOSH_TEST_DEPS) $(ACOSH_B)/hwy_acosh.o $(ACOSH_EXTRAS) common/include/hwymath/accuracy.h
	$(CXX) $(CXXFLAGS) $(INC) $< $(ACOSH_B)/hwy_acosh.o $(ACOSH_OBJS) $(ACOSH_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) -lm -o $@

# --- entry points ---
acosh-all: $(ACOSH_PROGS)
acosh-test: $(ACOSH_PROGS)
	@set -e; for t in $(ACOSH_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ACOSH_B)/hwy_test 2000000 --no-bench
acosh-hwy: $(ACOSH_B)/hwy_test
	./$(ACOSH_B)/hwy_test
acosh-hashes: $(ACOSH_B)/hwy_test
	$(RUN) ./$(ACOSH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 AVX2, bit for bit (pixi run -e cross make acosh-xcheck)
acosh-xcheck: $(ACOSH_GEN)
	functions/asinh/tools/xcheck.sh acosh $(PREFIX) $(PREFIX)
ALL_TARGETS    += acosh-all
TEST_TARGETS   += acosh-test
HWY_TARGETS    += acosh-hwy
HASHES_TARGETS += acosh-hashes
XCHECK_TARGETS += acosh-xcheck

ifeq ($(ARCH),x86_64)
acosh-insns: $(ACOSH_B)/svml_acosh_d_ha.o $(ACOSH_B)/svml_asm.o
	common/tools/compare_insns.sh $(ACOSH_B)/svml_acosh_d_ha.o svml_acosh8_ha $(ACOSH_B)/svml_asm.o __svml_acosh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly and of every Highway target (AVX-512 ones
# included) vs the translation
acosh-avx512-test: $(ACOSH_B)/test_acosh8 $(ACOSH_B)/hwy_test
	$(RUN) ./$(ACOSH_B)/test_acosh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ACOSH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += acosh-insns
AVX512_TARGETS += acosh-avx512-test
endif

.PHONY: acosh-all acosh-test acosh-hwy acosh-hashes acosh-xcheck acosh-insns acosh-avx512-test
