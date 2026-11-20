# log10, double precision: bit-identical C/AVX-512 translation of svml_z0_log10_d_ha.s and a Google Highway version
# that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths are relative
# to the repository root. Targets: log10-all log10-test log10-hwy log10-hashes log10-xcheck, and on x86-64 also log10-insns
# log10-avx512-test. Shares the constant generator, the Highway emulation helpers, the test driver and the
# comparison implementations with log (functions/log; see its function.mk).

LOG10     := functions/log10
LOG10_B   := build/log10
LOG10_ASM := $(SVML_ASM)/svml_z0_log10_d_ha.s
LOG10_H   := $(LOG10)/log10_d.h
# shared by the log family
LOG10F_GEN     := functions/log/tools/gen_svml_constants.py common/tools/extract_tables.py
LOG10F_HWY_H   := functions/log/hwy/log_hwy-inl.h functions/log/hwy/log_tables.h functions/log/hwy/rcp14_thresholds.h
LOG10F_TEST_H  := functions/log/tests/log_family_test.h functions/log/tests/extras.h common/include/hwymath/accuracy.h \
                common/include/hwymath/test_util.h

# --- generated constants (committed; regenerated when the generator changes) ---
LOG10_GEN := $(LOG10)/svml/svml_log10_d_ha_constants.inc $(LOG10)/hwy/hwy_log10_constants.inc
$(LOG10)/svml/svml_log10_d_ha_constants.inc: $(LOG10F_GEN) $(LOG10_ASM)
	python3 -I $< log10 $(LOG10_ASM) > $@
$(LOG10)/hwy/hwy_log10_constants.inc: $(LOG10F_GEN) $(LOG10_ASM)
	python3 -I $< log10 $(LOG10_ASM) --cxx > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
LOG10_OBJS   := $(LOG10_B)/svml_log10_d_ha.o $(LOG10_B)/svml_asm.o
LOG10_EXTRAS := $(LOG10_B)/extras_x86.o
LOG10_TESTS  := $(LOG10_B)/test_rare $(LOG10_B)/test_log108
else
LOG10_OBJS   := $(LOG10_B)/svml_log10_d_ha.o
LOG10_EXTRAS := $(LOG10_B)/extras_arm.o $(LOG10_B)/extras_arm_sve.o
LOG10_TESTS  :=
endif
LOG10_PROGS := $(LOG10_TESTS) $(LOG10_B)/hwy_test

$(LOG10_B):
	mkdir -p $@

# --- C translation ---
$(LOG10_B)/svml_log10_d_ha.o: $(LOG10)/svml/svml_log10_d_ha.c $(LOG10_H) $(LOG10_GEN) | $(LOG10_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(LOG10_B)/svml_asm.o: $(LOG10_ASM) | $(LOG10_B)
	sed 's/^__svml_dlog10_ha_cout_rare_internal:/\t.globl __svml_dlog10_ha_cout_rare_internal\n&/' $< > $(LOG10_B)/svml_asm.s
	$(CC) -c $(LOG10_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(LOG10_B)/extras_x86.o: functions/log/tests/extras_x86.c functions/log/tests/extras.h | $(LOG10_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -DFN=log10 -c $< -o $@
$(LOG10_B)/extras_arm.o: functions/log/tests/extras_arm.c functions/log/tests/extras.h | $(LOG10_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log10 -c $< -o $@
$(LOG10_B)/extras_arm_sve.o: functions/log/tests/extras_arm_sve.c | $(LOG10_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log10 -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
LOG10_TEST_DEPS := $(LOG10)/tests/inputs.h common/include/hwymath/test_util.h $(LOG10_OBJS)
$(LOG10_B)/test_%: $(LOG10)/tests/test_%.c $(LOG10_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(LOG10_OBJS) -lm -o $@

# --- Highway version ---
$(LOG10_B)/hwy_log10.o: $(LOG10)/hwy/hwy_log10.cc $(LOG10)/hwy/hwy_log10_constants.h $(LOG10F_HWY_H) $(LOG10_GEN) | $(LOG10_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(LOG10_B)/hwy_test: $(LOG10)/tests/hwy_test.cc $(LOG10F_TEST_H) $(LOG10_TEST_DEPS) $(LOG10_B)/hwy_log10.o $(LOG10_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LOG10_B)/hwy_log10.o $(LOG10_OBJS) $(LOG10_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
log10-all: $(LOG10_PROGS)
log10-test: $(LOG10_TESTS) $(LOG10_B)/hwy_test
	@set -e; for t in $(LOG10_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(LOG10_B)/hwy_test 2000000 --no-bench
log10-hwy: $(LOG10_B)/hwy_test
	./$(LOG10_B)/hwy_test
log10-hashes: $(LOG10_B)/hwy_test
	$(RUN) ./$(LOG10_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# Highway on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (pixi run -e cross make log10-xcheck)
log10-xcheck: $(LOG10_GEN)
	./functions/log/tools/xcheck.sh log10 $(PREFIX) $(PREFIX)
ALL_TARGETS    += log10-all
TEST_TARGETS   += log10-test
HWY_TARGETS    += log10-hwy
HASHES_TARGETS += log10-hashes
XCHECK_TARGETS += log10-xcheck

ifeq ($(ARCH),x86_64)
log10-insns: $(LOG10_B)/svml_log10_d_ha.o $(LOG10_B)/svml_asm.o
	common/tools/compare_insns.sh $(LOG10_B)/svml_log10_d_ha.o svml_log108_ha $(LOG10_B)/svml_asm.o __svml_log108_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly; Highway targets vs the translation
log10-avx512-test: $(LOG10_B)/test_log108 $(LOG10_B)/hwy_test
	$(RUN) ./$(LOG10_B)/test_log108 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(LOG10_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS   += log10-insns
AVX512_TARGETS  += log10-avx512-test
endif

.PHONY: log10-all log10-test log10-hwy log10-hashes log10-xcheck log10-insns log10-avx512-test
