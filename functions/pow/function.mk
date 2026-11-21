# pow, double precision: bit-identical C/AVX-512 translation of svml_z0_pow_d_ha.s and a Google Highway
# version of the same algorithm. Included by the top-level Makefile; paths are relative to the repository
# root. Targets: pow-all pow-test pow-hwy pow-hashes pow-xcheck, and on x86-64 also pow-insns pow-avx512-test
# pow-rcp14-table.

POW     := functions/pow
POW_B   := build/pow
POW_ASM := $(SVML_ASM)/svml_z0_pow_d_ha.s
POW_H   := $(POW)/pow_d.h $(POW)/svml/svml_pow_d_ha_constants.inc $(POW)/hwy/rcp14_switch.inc

# --- generated constants (committed; regenerated when their generator changes) ---
$(POW)/svml/svml_pow_d_ha_constants.inc: $(POW)/tools/gen_svml_constants.py common/tools/extract_tables.py $(POW_ASM)
	python3 -I $< $(POW_ASM) > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
POW_OBJS   := $(POW_B)/svml_pow_d_ha.o $(POW_B)/svml_asm.o
POW_EXTRAS := $(POW_B)/extras_x86.o
POW_TESTS  := $(POW_B)/test_rare $(POW_B)/test_pow8
else
POW_OBJS   := $(POW_B)/svml_pow_d_ha.o
POW_EXTRAS := $(POW_B)/extras_arm.o
POW_TESTS  :=
endif
POW_PROGS := $(POW_TESTS) $(POW_B)/hwy_test

$(POW_B):
	mkdir -p $@

$(POW_B)/svml_pow_d_ha.o: $(POW)/svml/svml_pow_d_ha.c $(POW_H) | $(POW_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(POW_B)/svml_asm.o: $(POW_ASM) | $(POW_B)
	sed 's/^__svml_dpow_ha_cout_rare_internal:/\t.globl __svml_dpow_ha_cout_rare_internal\n&/' $< > $(POW_B)/svml_asm.s
	$(CC) -c $(POW_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(POW_B)/extras_x86.o: $(POW)/tests/extras_x86.c $(POW)/tests/extras.h $(POW_H) | $(POW_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(POW_B)/extras_arm.o: $(POW)/tests/extras_arm.c $(POW)/tests/extras.h | $(POW_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@

# --- tests of the C translation (x86-64) ---
POW_TEST_DEPS := $(POW)/tests/inputs.h common/include/hwymath/test_util.h $(POW_OBJS)
$(POW_B)/test_%: $(POW)/tests/test_%.c $(POW_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(POW_OBJS) -lm -o $@

# switch points of the rounded vrcp14pd (needs AVX-512: run with RUN="$SDE -spr --"); committed output
$(POW_B)/gen_rcp14_table: $(POW)/tools/gen_rcp14_table.c | $(POW_B)
	$(CC) $(CFLAGS) $< -o $@
pow-rcp14-table: $(POW_B)/gen_rcp14_table
	$(RUN) ./$< > $(POW_B)/rcp14_switch.inc
	diff $(POW_B)/rcp14_switch.inc $(POW)/hwy/rcp14_switch.inc && echo "rcp14 switch table: identical"

# --- Highway version ---
$(POW_B)/hwy_pow.o: $(POW)/hwy/hwy_pow.cc $(POW)/hwy/hwy_pow_constants.h $(POW_H) | $(POW_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(POW_B)/hwy_test: $(POW)/tests/hwy_test.cc $(POW)/tests/extras.h common/include/hwymath/accuracy.h $(POW_TEST_DEPS) \
                   $(POW_B)/hwy_pow.o $(POW_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(POW_B)/hwy_pow.o $(POW_OBJS) $(POW_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
pow-all: $(POW_PROGS)
pow-test: $(POW_TESTS) $(POW_B)/hwy_test
	@set -e; for t in $(POW_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(POW_B)/hwy_test 2000000 --no-bench
pow-hwy: $(POW_B)/hwy_test
	./$(POW_B)/hwy_test
pow-hashes: $(POW_B)/hwy_test
	$(RUN) ./$(POW_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
pow-xcheck: $(POW)/svml/svml_pow_d_ha_constants.inc
	XC_FIXED="sve_fixed256:-march=armv8.2-a+sve -msve-vector-bits=256:32;sve2_fixed128:-march=armv9-a+sve2 -msve-vector-bits=128:16" \
	  common/tools/xcheck.sh pow $(PREFIX) $(PREFIX) 1000000
ALL_TARGETS    += pow-all
TEST_TARGETS   += pow-test
HWY_TARGETS    += pow-hwy
HASHES_TARGETS += pow-hashes
XCHECK_TARGETS += pow-xcheck

ifeq ($(ARCH),x86_64)
pow-insns: $(POW_B)/svml_pow_d_ha.o $(POW_B)/svml_asm.o
	common/tools/compare_insns.sh $(POW_B)/svml_pow_d_ha.o svml_pow8_ha $(POW_B)/svml_asm.o __svml_pow8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
pow-avx512-test: $(POW_B)/test_pow8 $(POW_B)/hwy_test
	$(RUN) ./$(POW_B)/test_pow8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(POW_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS   += pow-insns
AVX512_TARGETS  += pow-avx512-test
endif

.PHONY: pow-all pow-test pow-hwy pow-hashes pow-xcheck pow-insns pow-avx512-test pow-rcp14-table
