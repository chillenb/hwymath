# log2, double precision: bit-identical C/AVX-512 translation of svml_z0_log2_d_ha.s and a Google Highway version
# that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths are relative
# to the repository root. Targets: log2-all log2-test log2-hwy log2-hashes log2-xcheck, and on x86-64 also log2-insns
# log2-avx512-test. Shares the constant generator, the Highway emulation helpers, the test driver and the
# comparison implementations with log (functions/log; see its function.mk).

LOG2     := functions/log2
LOG2_B   := build/log2
LOG2_ASM := $(SVML_ASM)/svml_z0_log2_d_ha.s
LOG2_H   := $(LOG2)/log2_d.h
# shared by the log family
LOG2F_GEN     := functions/log/tools/gen_svml_constants.py common/tools/extract_tables.py
LOG2F_HWY_H   := functions/log/hwy/log_hwy-inl.h functions/log/hwy/log_tables.h functions/log/hwy/rcp14_thresholds.h
LOG2F_TEST_H  := functions/log/tests/log_family_test.h functions/log/tests/extras.h common/include/hwymath/accuracy.h \
                common/include/hwymath/test_util.h

# --- generated constants (committed; regenerated when the generator changes) ---
LOG2_GEN := $(LOG2)/svml/svml_log2_d_ha_constants.inc $(LOG2)/hwy/hwy_log2_constants.inc
$(LOG2)/svml/svml_log2_d_ha_constants.inc: $(LOG2F_GEN) $(LOG2_ASM)
	python3 -I $< log2 $(LOG2_ASM) > $@
$(LOG2)/hwy/hwy_log2_constants.inc: $(LOG2F_GEN) $(LOG2_ASM)
	python3 -I $< log2 $(LOG2_ASM) --cxx > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
LOG2_OBJS   := $(LOG2_B)/svml_log2_d_ha.o $(LOG2_B)/svml_asm.o
LOG2_EXTRAS := $(LOG2_B)/extras_x86.o
LOG2_TESTS  := $(LOG2_B)/test_rare $(LOG2_B)/test_log28
else
LOG2_OBJS   := $(LOG2_B)/svml_log2_d_ha.o
LOG2_EXTRAS := $(LOG2_B)/extras_arm.o $(LOG2_B)/extras_arm_sve.o
LOG2_TESTS  :=
endif
LOG2_PROGS := $(LOG2_TESTS) $(LOG2_B)/hwy_test

$(LOG2_B):
	mkdir -p $@

# --- C translation ---
$(LOG2_B)/svml_log2_d_ha.o: $(LOG2)/svml/svml_log2_d_ha.c $(LOG2_H) $(LOG2_GEN) | $(LOG2_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(LOG2_B)/svml_asm.o: $(LOG2_ASM) | $(LOG2_B)
	sed 's/^__svml_dlog2_ha_cout_rare_internal:/\t.globl __svml_dlog2_ha_cout_rare_internal\n&/' $< > $(LOG2_B)/svml_asm.s
	$(CC) -c $(LOG2_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(LOG2_B)/extras_x86.o: functions/log/tests/extras_x86.c functions/log/tests/extras.h | $(LOG2_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -DFN=log2 -c $< -o $@
$(LOG2_B)/extras_arm.o: functions/log/tests/extras_arm.c functions/log/tests/extras.h | $(LOG2_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log2 -c $< -o $@
$(LOG2_B)/extras_arm_sve.o: functions/log/tests/extras_arm_sve.c | $(LOG2_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log2 -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
LOG2_TEST_DEPS := $(LOG2)/tests/inputs.h common/include/hwymath/test_util.h $(LOG2_OBJS)
$(LOG2_B)/test_%: $(LOG2)/tests/test_%.c $(LOG2_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(LOG2_OBJS) -lm -o $@

# --- Highway version ---
$(LOG2_B)/hwy_log2.o: $(LOG2)/hwy/hwy_log2.cc $(LOG2)/hwy/hwy_log2_constants.h $(LOG2F_HWY_H) $(LOG2_GEN) | $(LOG2_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(LOG2_B)/hwy_test: $(LOG2)/tests/hwy_test.cc $(LOG2F_TEST_H) $(LOG2_TEST_DEPS) $(LOG2_B)/hwy_log2.o $(LOG2_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LOG2_B)/hwy_log2.o $(LOG2_OBJS) $(LOG2_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
log2-all: $(LOG2_PROGS)
log2-test: $(LOG2_TESTS) $(LOG2_B)/hwy_test
	@set -e; for t in $(LOG2_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(LOG2_B)/hwy_test 2000000 --no-bench
log2-hwy: $(LOG2_B)/hwy_test
	./$(LOG2_B)/hwy_test
log2-hashes: $(LOG2_B)/hwy_test
	$(RUN) ./$(LOG2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# Highway on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (pixi run -e cross make log2-xcheck)
log2-xcheck: $(LOG2_GEN)
	common/tools/xcheck.sh log2 $(PREFIX) $(PREFIX)
ALL_TARGETS    += log2-all
TEST_TARGETS   += log2-test
HWY_TARGETS    += log2-hwy
HASHES_TARGETS += log2-hashes
XCHECK_TARGETS += log2-xcheck

ifeq ($(ARCH),x86_64)
log2-insns: $(LOG2_B)/svml_log2_d_ha.o $(LOG2_B)/svml_asm.o
	common/tools/compare_insns.sh $(LOG2_B)/svml_log2_d_ha.o svml_log28_ha $(LOG2_B)/svml_asm.o __svml_log28_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly; Highway targets vs the translation
log2-avx512-test: $(LOG2_B)/test_log28 $(LOG2_B)/hwy_test
	$(RUN) ./$(LOG2_B)/test_log28 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(LOG2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS   += log2-insns
AVX512_TARGETS  += log2-avx512-test
endif

.PHONY: log2-all log2-test log2-hwy log2-hashes log2-xcheck log2-insns log2-avx512-test
