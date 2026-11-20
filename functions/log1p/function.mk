# log1p, double precision: bit-identical C/AVX-512 translation of svml_z0_log1p_d_ha.s and a Google Highway version
# that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths are relative
# to the repository root. Targets: log1p-all log1p-test log1p-hwy log1p-hashes log1p-xcheck, and on x86-64 also log1p-insns
# log1p-avx512-test. Shares the constant generator, the Highway emulation helpers, the test driver and the
# comparison implementations with log (functions/log; see its function.mk).

LOG1P     := functions/log1p
LOG1P_B   := build/log1p
LOG1P_ASM := $(SVML_ASM)/svml_z0_log1p_d_ha.s
LOG1P_H   := $(LOG1P)/log1p_d.h
# shared by the log family
LOG1PF_GEN     := functions/log/tools/gen_svml_constants.py common/tools/extract_tables.py
LOG1PF_HWY_H   := functions/log/hwy/log_hwy-inl.h functions/log/hwy/log_tables.h functions/log/hwy/rcp14_thresholds.h
LOG1PF_TEST_H  := functions/log/tests/log_family_test.h functions/log/tests/extras.h common/include/hwymath/accuracy.h \
                common/include/hwymath/test_util.h

# --- generated constants (committed; regenerated when the generator changes) ---
LOG1P_GEN := $(LOG1P)/svml/svml_log1p_d_ha_constants.inc $(LOG1P)/hwy/hwy_log1p_constants.inc
$(LOG1P)/svml/svml_log1p_d_ha_constants.inc: $(LOG1PF_GEN) $(LOG1P_ASM)
	python3 -I $< log1p $(LOG1P_ASM) > $@
$(LOG1P)/hwy/hwy_log1p_constants.inc: $(LOG1PF_GEN) $(LOG1P_ASM)
	python3 -I $< log1p $(LOG1P_ASM) --cxx > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
LOG1P_OBJS   := $(LOG1P_B)/svml_log1p_d_ha.o $(LOG1P_B)/svml_asm.o
LOG1P_EXTRAS := $(LOG1P_B)/extras_x86.o
LOG1P_TESTS  := $(LOG1P_B)/test_rare $(LOG1P_B)/test_log1p8
else
LOG1P_OBJS   := $(LOG1P_B)/svml_log1p_d_ha.o
LOG1P_EXTRAS := $(LOG1P_B)/extras_arm.o $(LOG1P_B)/extras_arm_sve.o
LOG1P_TESTS  :=
endif
LOG1P_PROGS := $(LOG1P_TESTS) $(LOG1P_B)/hwy_test

$(LOG1P_B):
	mkdir -p $@

# --- C translation ---
$(LOG1P_B)/svml_log1p_d_ha.o: $(LOG1P)/svml/svml_log1p_d_ha.c $(LOG1P_H) $(LOG1P_GEN) | $(LOG1P_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(LOG1P_B)/svml_asm.o: $(LOG1P_ASM) | $(LOG1P_B)
	sed 's/^__svml_dlog1p_ha_cout_rare_internal:/\t.globl __svml_dlog1p_ha_cout_rare_internal\n&/' $< > $(LOG1P_B)/svml_asm.s
	$(CC) -c $(LOG1P_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(LOG1P_B)/extras_x86.o: functions/log/tests/extras_x86.c functions/log/tests/extras.h | $(LOG1P_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -DFN=log1p -c $< -o $@
$(LOG1P_B)/extras_arm.o: functions/log/tests/extras_arm.c functions/log/tests/extras.h | $(LOG1P_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log1p -c $< -o $@
$(LOG1P_B)/extras_arm_sve.o: functions/log/tests/extras_arm_sve.c | $(LOG1P_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log1p -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
LOG1P_TEST_DEPS := $(LOG1P)/tests/inputs.h common/include/hwymath/test_util.h $(LOG1P_OBJS)
$(LOG1P_B)/test_%: $(LOG1P)/tests/test_%.c $(LOG1P_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(LOG1P_OBJS) -lm -o $@

# --- Highway version ---
$(LOG1P_B)/hwy_log1p.o: $(LOG1P)/hwy/hwy_log1p.cc $(LOG1P)/hwy/hwy_log1p_constants.h $(LOG1PF_HWY_H) $(LOG1P_GEN) | $(LOG1P_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(LOG1P_B)/hwy_test: $(LOG1P)/tests/hwy_test.cc $(LOG1PF_TEST_H) $(LOG1P_TEST_DEPS) $(LOG1P_B)/hwy_log1p.o $(LOG1P_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LOG1P_B)/hwy_log1p.o $(LOG1P_OBJS) $(LOG1P_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
log1p-all: $(LOG1P_PROGS)
log1p-test: $(LOG1P_TESTS) $(LOG1P_B)/hwy_test
	@set -e; for t in $(LOG1P_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(LOG1P_B)/hwy_test 2000000 --no-bench
log1p-hwy: $(LOG1P_B)/hwy_test
	./$(LOG1P_B)/hwy_test
log1p-hashes: $(LOG1P_B)/hwy_test
	$(RUN) ./$(LOG1P_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# Highway on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (pixi run -e cross make log1p-xcheck)
log1p-xcheck: $(LOG1P_GEN)
	./functions/log/tools/xcheck.sh log1p $(PREFIX) $(PREFIX)
ALL_TARGETS    += log1p-all
TEST_TARGETS   += log1p-test
HWY_TARGETS    += log1p-hwy
HASHES_TARGETS += log1p-hashes
XCHECK_TARGETS += log1p-xcheck

ifeq ($(ARCH),x86_64)
log1p-insns: $(LOG1P_B)/svml_log1p_d_ha.o $(LOG1P_B)/svml_asm.o
	common/tools/compare_insns.sh $(LOG1P_B)/svml_log1p_d_ha.o svml_log1p8_ha $(LOG1P_B)/svml_asm.o __svml_log1p8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly; Highway targets vs the translation
log1p-avx512-test: $(LOG1P_B)/test_log1p8 $(LOG1P_B)/hwy_test
	$(RUN) ./$(LOG1P_B)/test_log1p8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(LOG1P_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS   += log1p-insns
AVX512_TARGETS  += log1p-avx512-test
endif

.PHONY: log1p-all log1p-test log1p-hwy log1p-hashes log1p-xcheck log1p-insns log1p-avx512-test
