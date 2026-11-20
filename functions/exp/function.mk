# exp, double precision: bit-identical C/AVX-512 translation of svml_z0_exp_d_ha.s, AVX2 ports and
# variants, and a Google Highway version. Included by the top-level Makefile; paths are relative to the
# repository root. Targets: exp-all exp-test exp-hwy exp-hashes exp-xcheck, and on x86-64 also
# exp-compare exp-bench exp-insns exp-avx512-test.

EXP     := functions/exp
EXP_B   := build/exp
EXP_ASM := $(SVML_ASM)/svml_z0_exp_d_ha.s
EXP_H   := $(EXP)/exp_d.h

# --- generated constants (committed; regenerated when their generators change) ---
EXP_FAST_K   ?= 3
EXP_FAST_DEG ?= 7
EXP_HP_K     ?= 2
EXP_HP_DEG   ?= 9
EXP_GEN := $(EXP)/svml/svml_exp_d_ha_constants.inc $(EXP)/avx2/exp_d_fast_constants.inc \
           $(EXP)/avx2/exp_d_hp_constants.inc $(EXP)/hwy/hwy_exp_constants.inc

$(EXP)/svml/svml_exp_d_ha_constants.inc: $(EXP)/tools/gen_svml_constants.py common/tools/extract_tables.py $(EXP_ASM)
	python3 -I $< $(EXP_ASM) > $@
$(EXP)/avx2/exp_d_fast_constants.inc: $(EXP)/tools/fit_poly.py
	python3 -I $< $(EXP_FAST_K) $(EXP_FAST_DEG) > $@
$(EXP)/avx2/exp_d_hp_constants.inc: $(EXP)/tools/fit_poly.py
	python3 -I $< $(EXP_HP_K) $(EXP_HP_DEG) hp > $@
$(EXP)/hwy/hwy_exp_constants.inc: $(EXP)/tools/fit_poly.py
	( python3 -I $< 2 9 hp | sed -E 's/\bh(Th|Tlr|C[0-9]+|Shifter|ShifterBase|K|Ln2Hi|Ln2Lo)\b/h2\1/g'; \
	  python3 -I $< 3 8 hp | sed -E 's/\bh(Th|Tlr|C[0-9]+|Shifter|ShifterBase|K|Ln2Hi|Ln2Lo)\b/h3\1/g' ) \
	  | sed 's/^static const double \(h[23]T\)/static constexpr double \1/' > $@

# --- what gets built on this architecture ---
ifeq ($(ARCH),x86_64)
EXP_OBJS   := $(EXP_B)/svml_exp_d_ha.o $(EXP_B)/svml_exp_d_ha_avx2.o $(EXP_B)/exp_d_avx2_fast.o \
              $(EXP_B)/exp_d_avx2_hp.o $(EXP_B)/svml_asm.o
EXP_EXTRAS := $(EXP_B)/extras_x86.o
EXP_TESTS  := $(EXP_B)/test_rare $(EXP_B)/test_model $(EXP_B)/test_exp8
EXP_PROGS  := $(EXP_TESTS) $(EXP_B)/hwy_test $(if $(INTEL_LIB),$(EXP_B)/compare_avx2 $(EXP_B)/bench)
else
EXP_OBJS   := $(EXP_B)/svml_exp_d_ha.o
EXP_EXTRAS := $(EXP_B)/extras_arm.o $(EXP_B)/extras_arm_sve.o
EXP_TESTS  :=
EXP_PROGS  := $(EXP_B)/hwy_test
endif

$(EXP_B):
	mkdir -p $@

# --- C kernels ---
$(EXP_B)/%.o: $(EXP)/svml/%.c $(EXP_H) $(EXP_GEN) | $(EXP_B)
	$(CC) $(CFLAGS) -c $< -o $@
$(EXP_B)/%.o: $(EXP)/avx2/%.c $(EXP_H) $(EXP_GEN) | $(EXP_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(EXP_B)/svml_asm.o: $(EXP_ASM) | $(EXP_B)
	sed 's/^__svml_dexp_ha_cout_rare_internal:/\t.globl __svml_dexp_ha_cout_rare_internal\n&/' $< > $(EXP_B)/svml_asm.s
	$(CC) -c $(EXP_B)/svml_asm.s -o $@

# --- comparison implementations ---
$(EXP_B)/extras_x86.o: $(EXP)/tests/extras_x86.c $(EXP)/tests/extras.h $(EXP_H) | $(EXP_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(EXP_B)/extras_arm.o: $(EXP)/tests/extras_arm.c $(EXP)/tests/extras.h | $(EXP_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
$(EXP_B)/extras_arm_sve.o: $(EXP)/tests/extras_arm_sve.c | $(EXP_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests and benchmarks of the C kernels (x86-64) ---
EXP_TEST_DEPS := $(EXP)/tests/inputs.h common/include/hwymath/test_util.h $(EXP_OBJS)
$(EXP_B)/test_model: $(EXP)/tests/test_model.c $(EXP_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) -frounding-math $< $(EXP_OBJS) -lm -o $@
$(EXP_B)/compare_avx2 $(EXP_B)/bench: $(EXP_B)/%: $(EXP)/tests/%.c $(EXP_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $< $(EXP_OBJS) $(LIBPATH) $(SLEEF_LIB) $(INTEL_LIB) -lmvec -lm -o $@
$(EXP_B)/test_%: $(EXP)/tests/test_%.c $(EXP_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(EXP_OBJS) -lm -o $@

# --- Highway version ---
$(EXP_B)/hwy_exp.o: $(EXP)/hwy/hwy_exp.cc $(EXP)/hwy/hwy_exp_constants.h $(EXP_GEN) | $(EXP_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(EXP_B)/hwy_test: $(EXP)/tests/hwy_test.cc $(EXP)/tests/extras.h $(EXP_TEST_DEPS) $(EXP_B)/hwy_exp.o $(EXP_EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< $(EXP_B)/hwy_exp.o $(EXP_OBJS) $(EXP_EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) \
	  $(INTEL_LIB) $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
exp-all: $(EXP_PROGS)
exp-test: $(EXP_TESTS) $(EXP_B)/hwy_test
	@set -e; for t in $(EXP_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(EXP_B)/hwy_test 2000000 --no-bench
exp-hwy: $(EXP_B)/hwy_test
	./$(EXP_B)/hwy_test
exp-hashes: $(EXP_B)/hwy_test
	$(RUN) ./$(EXP_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
exp-xcheck: $(EXP_GEN)
	./$(EXP)/tools/xcheck.sh $(PREFIX) $(PREFIX)
ALL_TARGETS    += exp-all
TEST_TARGETS   += exp-test
HWY_TARGETS    += exp-hwy
HASHES_TARGETS += exp-hashes
XCHECK_TARGETS += exp-xcheck

ifeq ($(ARCH),x86_64)
exp-compare: $(EXP_B)/compare_avx2
	./$(EXP_B)/compare_avx2
exp-bench: $(EXP_B)/bench
	./$(EXP_B)/bench
exp-insns: $(EXP_B)/svml_exp_d_ha.o $(EXP_B)/svml_asm.o
	common/tools/compare_insns.sh $(EXP_B)/svml_exp_d_ha.o svml_exp8_ha $(EXP_B)/svml_asm.o __svml_exp8_ha
# AVX-512: bit-exactness of the translation (and AVX2 port) vs Intel's assembly, Highway AVX-512 targets
exp-avx512-test: $(EXP_B)/test_exp8 $(EXP_B)/hwy_test
	$(RUN) ./$(EXP_B)/test_exp8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(EXP_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
COMPARE_TARGETS += exp-compare
BENCH_TARGETS   += exp-bench
INSNS_TARGETS   += exp-insns
AVX512_TARGETS  += exp-avx512-test
endif

.PHONY: exp-all exp-test exp-hwy exp-hashes exp-xcheck exp-compare exp-bench exp-insns exp-avx512-test
