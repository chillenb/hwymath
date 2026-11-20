# log, double precision: bit-identical C/AVX-512 translation of svml_z0_log_d_ha.s and a Google Highway version
# that is bit-identical to it on every target with FMA. Included by the top-level Makefile; paths are relative
# to the repository root. Targets: log-all log-test log-hwy log-hashes log-xcheck, and on x86-64 also log-insns
# log-avx512-test log-rcp14-check.
# functions/log also holds the parts shared with log2, log10 and log1p: tools/gen_svml_constants.py,
# hwy/log_hwy-inl.h (emulation of the AVX-512 operations), hwy/rcp14_thresholds.h, tools/verify_rcp14.c,
# tests/log_family_test.h and tests/extras*.

LOG     := functions/log
LOG_B   := build/log
LOG_ASM := $(SVML_ASM)/svml_z0_log_d_ha.s
LOG_H   := $(LOG)/log_d.h
# shared by the log family
LOGF_GEN     := $(LOG)/tools/gen_svml_constants.py common/tools/extract_tables.py
LOGF_HWY_H   := $(LOG)/hwy/log_hwy-inl.h $(LOG)/hwy/log_tables.h $(LOG)/hwy/rcp14_thresholds.h
LOGF_TEST_H  := $(LOG)/tests/log_family_test.h $(LOG)/tests/extras.h common/include/hwymath/accuracy.h \
                common/include/hwymath/test_util.h

# --- generated constants (committed; regenerated when the generator changes) ---
LOG_GEN := $(LOG)/svml/svml_log_d_ha_constants.inc $(LOG)/hwy/hwy_log_constants.inc
$(LOG)/svml/svml_log_d_ha_constants.inc: $(LOGF_GEN) $(LOG_ASM)
	python3 -I $< log $(LOG_ASM) > $@
$(LOG)/hwy/hwy_log_constants.inc: $(LOGF_GEN) $(LOG_ASM)
	python3 -I $< log $(LOG_ASM) --cxx > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
LOG_OBJS   := $(LOG_B)/svml_log_d_ha.o $(LOG_B)/svml_asm.o
LOG_EXTRAS := $(LOG_B)/extras_x86.o
LOG_TESTS  := $(LOG_B)/test_rare $(LOG_B)/test_log8 $(LOG_B)/verify_rcp14
else
LOG_OBJS   := $(LOG_B)/svml_log_d_ha.o
LOG_EXTRAS := $(LOG_B)/extras_arm.o $(LOG_B)/extras_arm_sve.o
LOG_TESTS  :=
endif
LOG_PROGS := $(LOG_TESTS) $(LOG_B)/hwy_test

$(LOG_B):
	mkdir -p $@

# --- C translation ---
$(LOG_B)/svml_log_d_ha.o: $(LOG)/svml/svml_log_d_ha.c $(LOG_H) $(LOG_GEN) | $(LOG_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(LOG_B)/svml_asm.o: $(LOG_ASM) | $(LOG_B)
	sed 's/^__svml_dlog_ha_cout_rare_internal:/\t.globl __svml_dlog_ha_cout_rare_internal\n&/' $< > $(LOG_B)/svml_asm.s
	$(CC) -c $(LOG_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(LOG_B)/extras_x86.o: $(LOG)/tests/extras_x86.c $(LOG)/tests/extras.h | $(LOG_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -DFN=log -c $< -o $@
$(LOG_B)/extras_arm.o: $(LOG)/tests/extras_arm.c $(LOG)/tests/extras.h | $(LOG_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log -c $< -o $@
$(LOG_B)/extras_arm_sve.o: $(LOG)/tests/extras_arm_sve.c | $(LOG_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -DFN=log -march=armv8-a+sve -c $< -o $@

# --- tests of the C translation (x86-64) ---
LOG_TEST_DEPS := $(LOG)/tests/inputs.h common/include/hwymath/test_util.h $(LOG_OBJS)
$(LOG_B)/test_%: $(LOG)/tests/test_%.c $(LOG_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(LOG_OBJS) -lm -o $@
# the rounded vrcp14pd model (rcp14_thresholds.h) vs Intel's reference algorithm (rcp14_ref.h) and the CPU
$(LOG_B)/verify_rcp14: $(LOG)/tools/verify_rcp14.c $(LOG)/tools/rcp14_ref.h $(LOG)/hwy/rcp14_thresholds.h | $(LOG_B)
	$(CXX) $(CXXFLAGS) $(INC) -x c++ $< -lm -o $@

# --- Highway version ---
$(LOG_B)/hwy_log.o: $(LOG)/hwy/hwy_log.cc $(LOG)/hwy/hwy_log_constants.h $(LOGF_HWY_H) $(LOG_GEN) | $(LOG_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(LOG_B)/hwy_test: $(LOG)/tests/hwy_test.cc $(LOGF_TEST_H) $(LOG_TEST_DEPS) $(LOG_B)/hwy_log.o $(LOG_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(LOG_B)/hwy_log.o $(LOG_OBJS) $(LOG_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
log-all: $(LOG_PROGS)
log-test: $(LOG_TESTS) $(LOG_B)/hwy_test
	@set -e; for t in $(LOG_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(LOG_B)/hwy_test 2000000 --no-bench
log-hwy: $(LOG_B)/hwy_test
	./$(LOG_B)/hwy_test
log-hashes: $(LOG_B)/hwy_test
	$(RUN) ./$(LOG_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# Highway on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (pixi run -e cross make log-xcheck)
log-xcheck: $(LOG_GEN)
	./functions/log/tools/xcheck.sh log $(PREFIX) $(PREFIX)
ALL_TARGETS    += log-all
TEST_TARGETS   += log-test
HWY_TARGETS    += log-hwy
HASHES_TARGETS += log-hashes
XCHECK_TARGETS += log-xcheck

ifeq ($(ARCH),x86_64)
log-insns: $(LOG_B)/svml_log_d_ha.o $(LOG_B)/svml_asm.o
	common/tools/compare_insns.sh $(LOG_B)/svml_log_d_ha.o svml_log8_ha $(LOG_B)/svml_asm.o __svml_log8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly; Highway targets vs the translation
log-avx512-test: $(LOG_B)/test_log8 $(LOG_B)/hwy_test $(LOG_B)/verify_rcp14
	$(RUN) ./$(LOG_B)/test_log8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(LOG_B)/verify_rcp14 || [ $$? -eq 77 ]
	$(RUN) ./$(LOG_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
# the rounded-vrcp14pd model shared by the log family, against the CPU's vrcp14pd (skipped without AVX-512)
log-rcp14-check: $(LOG_B)/verify_rcp14
	$(RUN) ./$(LOG_B)/verify_rcp14 || [ $$? -eq 77 ]
INSNS_TARGETS   += log-insns
AVX512_TARGETS  += log-avx512-test
endif

.PHONY: log-all log-test log-hwy log-hashes log-xcheck log-insns log-avx512-test log-rcp14-check
