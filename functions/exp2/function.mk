# exp2, double precision: bit-identical C/AVX-512 translation of svml_z0_exp2_d_ha.s and a Google Highway
# version of the same algorithm (bit-identical on all FMA targets). Included by the top-level Makefile;
# paths are relative to the repository root. Targets: exp2-all exp2-test exp2-hwy exp2-hashes exp2-xcheck, and on
# x86-64 also exp2-insns exp2-avx512-test.

EXP2     := functions/exp2
EXP2_B   := build/exp2
EXP2_ASM := $(SVML_ASM)/svml_z0_exp2_d_ha.s
EXP2_H   := $(EXP2)/exp2_d.h
EXP2_GEN := $(EXP2)/svml/svml_exp2_d_ha_constants.inc

# generated constants (committed; regenerated when the generator changes)
$(EXP2_GEN): $(EXP2)/tools/gen_svml_constants.py common/tools/extract_tables.py $(EXP2_ASM)
	python3 -I $< $(EXP2_ASM) > $@

ifeq ($(ARCH),x86_64)
EXP2_OBJS  := $(EXP2_B)/svml_exp2_d_ha.o $(EXP2_B)/svml_asm.o
EXP2_TESTS := $(EXP2_B)/test_rare $(EXP2_B)/test_exp28
else
EXP2_OBJS  := $(EXP2_B)/svml_exp2_d_ha.o
EXP2_TESTS :=
endif

$(EXP2_B):
	mkdir -p $@

$(EXP2_B)/svml_exp2_d_ha.o: $(EXP2)/svml/svml_exp2_d_ha.c $(EXP2_H) $(EXP2_GEN) | $(EXP2_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(EXP2_B)/svml_asm.o: $(EXP2_ASM) | $(EXP2_B)
	sed 's/^__svml_dexp2_ha_cout_rare_internal:/\t.globl __svml_dexp2_ha_cout_rare_internal\n&/' $< > $(EXP2_B)/svml_asm.s
	$(CC) -c $(EXP2_B)/svml_asm.s -o $@

$(EXP2_B)/extras.o: $(EXP2)/tests/extras.c $(EXP2)/tests/extras.h $(EXP2_H) | $(EXP2_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@

EXP2_TEST_DEPS := $(EXP2)/tests/inputs.h common/include/hwymath/test_util.h $(EXP2_OBJS)
$(EXP2_B)/test_%: $(EXP2)/tests/test_%.c $(EXP2_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(EXP2_OBJS) -lm -o $@

# Highway version
$(EXP2_B)/hwy_exp2.o: $(EXP2)/hwy/hwy_exp2.cc $(EXP2)/hwy/hwy_exp2_constants.h $(EXP2_GEN) | $(EXP2_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(EXP2_B)/hwy_test: $(EXP2)/tests/hwy_test.cc $(EXP2)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(EXP2_TEST_DEPS) $(EXP2_B)/hwy_exp2.o $(EXP2_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(EXP2_B)/hwy_exp2.o $(EXP2_OBJS) $(EXP2_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

exp2-all: $(EXP2_TESTS) $(EXP2_B)/hwy_test
exp2-test: $(EXP2_TESTS) $(EXP2_B)/hwy_test
	@set -e; for t in $(EXP2_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(EXP2_B)/hwy_test 2000000 --no-bench
exp2-hwy: $(EXP2_B)/hwy_test
	./$(EXP2_B)/hwy_test
exp2-xcheck: $(EXP2_GEN)
	./$(EXP2)/tools/xcheck.sh $(PREFIX) $(PREFIX)
exp2-hashes: $(EXP2_B)/hwy_test
	$(RUN) ./$(EXP2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += exp2-all
TEST_TARGETS   += exp2-test
HWY_TARGETS    += exp2-hwy
HASHES_TARGETS += exp2-hashes
XCHECK_TARGETS += exp2-xcheck

ifeq ($(ARCH),x86_64)
exp2-insns: $(EXP2_B)/svml_exp2_d_ha.o $(EXP2_B)/svml_asm.o
	common/tools/compare_insns.sh $(EXP2_B)/svml_exp2_d_ha.o svml_exp28_ha $(EXP2_B)/svml_asm.o __svml_exp28_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway targets vs the translation
exp2-avx512-test: $(EXP2_B)/test_exp28 $(EXP2_B)/hwy_test
	$(RUN) ./$(EXP2_B)/test_exp28 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(EXP2_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += exp2-insns
AVX512_TARGETS += exp2-avx512-test
endif

.PHONY: exp2-all exp2-test exp2-hwy exp2-hashes exp2-xcheck exp2-insns exp2-avx512-test
