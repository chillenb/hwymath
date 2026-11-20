# sin, double precision: bit-identical C/AVX-512 translation of svml_z0_sin_d_ha.s and a Google Highway
# version. Also owns the large-argument reduction shared with cos and tan (svml/svml_trig_reduction.h,
# svml/svml_trig_reduction_data.inc, hwy/trig_hwy-inl.h), the sin/cos constant generator, the test
# templates (tests/*_tmpl.h, tests/inputs.h) and trig_common.mk. Included by the top-level Makefile;
# paths are relative to the repository root. Targets: sin-all sin-test sin-hwy sin-hashes sin-xcheck,
# and on x86-64 also sin-insns sin-avx512-test.

SIN     := functions/sin
SIN_B   := build/sin
SIN_ASM := $(SVML_ASM)/svml_z0_sin_d_ha.s
SIN_H   := $(SIN)/sin_d.h

include functions/sin/trig_common.mk

# --- generated constants (committed; regenerated when their generator changes) ---
SIN_GEN := $(SIN)/svml/svml_sin_d_ha_constants.inc $(TRIG_RED_INC)
$(SIN)/svml/svml_sin_d_ha_constants.inc: $(SIN)/tools/gen_svml_constants.py common/tools/extract_tables.py $(SIN_ASM)
	python3 -I $< sin $(SIN_ASM) > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
SIN_OBJS   := $(SIN_B)/svml_sin_d_ha.o $(SIN_B)/svml_asm.o
SIN_TESTS  := $(SIN_B)/test_rare $(SIN_B)/test_sin8
else
SIN_OBJS   := $(SIN_B)/svml_sin_d_ha.o
SIN_TESTS  :=
endif
SIN_PROGS := $(SIN_TESTS) $(SIN_B)/hwy_test

$(SIN_B):
	mkdir -p $@

$(SIN_B)/svml_sin_d_ha.o: $(SIN)/svml/svml_sin_d_ha.c $(SIN_H) $(SIN_GEN) $(TRIG_SVML_DEPS) | $(SIN_B)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(SIN_B)/svml_asm.o: $(SIN_ASM) | $(SIN_B)
	sed 's/^__svml_dsin_ha_cout_rare_internal:/\t.globl __svml_dsin_ha_cout_rare_internal\n&/' $< > $(SIN_B)/svml_asm.s
	$(CC) -c $(SIN_B)/svml_asm.s -o $@

# --- tests of the C translation (x86-64) ---
$(SIN_B)/test_%: $(SIN)/tests/test_%.c $(TRIG_TEST_DEPS) $(SIN_H) $(SIN_OBJS)
	$(CC) $(CFLAGS) $(INC) $< $(SIN_OBJS) -lm -o $@

# --- Highway version ---
$(SIN_B)/hwy_sin.o: $(SIN)/hwy/hwy_sin.cc $(SIN_GEN) $(TRIG_HWY_DEPS) | $(SIN_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(SIN_B)/extras.o: $(SIN)/tests/extras.c $(SIN)/tests/extras_tmpl.h $(SIN)/tests/extras_tmpl_decl.h $(SIN_H) | $(SIN_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(SIN_B)/hwy_test: $(SIN)/tests/hwy_test.cc $(TRIG_HWYTEST_DEPS) $(SIN_H) $(SIN_OBJS) $(SIN_B)/hwy_sin.o $(SIN_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(SIN_B)/hwy_sin.o $(SIN_OBJS) $(SIN_B)/extras.o $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
sin-all: $(SIN_PROGS)
sin-test: $(SIN_TESTS) $(SIN_B)/hwy_test
	@set -e; for t in $(SIN_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(SIN_B)/hwy_test 2000000 --no-bench
sin-hwy: $(SIN_B)/hwy_test
	./$(SIN_B)/hwy_test
sin-hashes: $(SIN_B)/hwy_test
	$(RUN) ./$(SIN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += sin-all
TEST_TARGETS   += sin-test
HWY_TARGETS    += sin-hwy
HASHES_TARGETS += sin-hashes
# Highway sin on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64 (run in the cross environment:
# pixi run -e cross make sin-xcheck)
sin-xcheck: $(SIN_GEN)
	./functions/sin/tools/xcheck.sh sin $(PREFIX) $(PREFIX)
XCHECK_TARGETS += sin-xcheck

ifeq ($(ARCH),x86_64)
sin-insns: $(SIN_B)/svml_sin_d_ha.o $(SIN_B)/svml_asm.o
	common/tools/compare_insns.sh $(SIN_B)/svml_sin_d_ha.o svml_sin8_ha $(SIN_B)/svml_asm.o __svml_sin8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
sin-avx512-test: $(SIN_B)/test_sin8 $(SIN_B)/hwy_test
	$(RUN) ./$(SIN_B)/test_sin8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(SIN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += sin-insns
AVX512_TARGETS += sin-avx512-test
endif

.PHONY: sin-all sin-test sin-hwy sin-hashes sin-xcheck sin-insns sin-avx512-test
