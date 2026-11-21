# expm1, double precision: bit-identical C/AVX-512 translation of svml_z0_expm1_d_ha.s and a Google Highway
# version of the same algorithm (bit-identical on all FMA targets). Included by the top-level Makefile;
# paths are relative to the repository root. Targets: expm1-all expm1-test expm1-hwy expm1-hashes expm1-xcheck, and on
# x86-64 also expm1-insns expm1-avx512-test.

EXPM1     := functions/expm1
EXPM1_B   := build/expm1
EXPM1_ASM := $(SVML_ASM)/svml_z0_expm1_d_ha.s
EXPM1_H   := $(EXPM1)/expm1_d.h
EXPM1_GEN := $(EXPM1)/svml/svml_expm1_d_ha_constants.inc

# generated constants (committed; regenerated when the generator changes)
$(EXPM1_GEN): $(EXPM1)/tools/gen_svml_constants.py common/tools/extract_tables.py $(EXPM1_ASM)
	python3 -I $< $(EXPM1_ASM) > $@

ifeq ($(ARCH),x86_64)
EXPM1_OBJS  := $(EXPM1_B)/svml_expm1_d_ha.o $(EXPM1_B)/svml_asm.o
EXPM1_TESTS := $(EXPM1_B)/test_rare $(EXPM1_B)/test_expm18
else
EXPM1_OBJS  := $(EXPM1_B)/svml_expm1_d_ha.o
EXPM1_TESTS :=
endif

$(EXPM1_B):
	mkdir -p $@

$(EXPM1_B)/svml_expm1_d_ha.o: $(EXPM1)/svml/svml_expm1_d_ha.c $(EXPM1_H) $(EXPM1_GEN) | $(EXPM1_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(EXPM1_B)/svml_asm.o: $(EXPM1_ASM) | $(EXPM1_B)
	sed 's/^__svml_dexpm1_ha_cout_rare_internal:/\t.globl __svml_dexpm1_ha_cout_rare_internal\n&/' $< > $(EXPM1_B)/svml_asm.s
	$(CC) -c $(EXPM1_B)/svml_asm.s -o $@

$(EXPM1_B)/extras.o: $(EXPM1)/tests/extras.c $(EXPM1)/tests/extras.h $(EXPM1_H) | $(EXPM1_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@

EXPM1_TEST_DEPS := $(EXPM1)/tests/inputs.h common/include/hwymath/test_util.h $(EXPM1_OBJS)
$(EXPM1_B)/test_%: $(EXPM1)/tests/test_%.c $(EXPM1_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(EXPM1_OBJS) -lm -o $@

# Highway version
$(EXPM1_B)/hwy_expm1.o: $(EXPM1)/hwy/hwy_expm1.cc $(EXPM1)/hwy/hwy_expm1_constants.h $(EXPM1_GEN) | $(EXPM1_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(EXPM1_B)/hwy_test: $(EXPM1)/tests/hwy_test.cc $(EXPM1)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(EXPM1_TEST_DEPS) $(EXPM1_B)/hwy_expm1.o $(EXPM1_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(EXPM1_B)/hwy_expm1.o $(EXPM1_OBJS) $(EXPM1_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

expm1-all: $(EXPM1_TESTS) $(EXPM1_B)/hwy_test
expm1-test: $(EXPM1_TESTS) $(EXPM1_B)/hwy_test
	@set -e; for t in $(EXPM1_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(EXPM1_B)/hwy_test 2000000 --no-bench
expm1-hwy: $(EXPM1_B)/hwy_test
	./$(EXPM1_B)/hwy_test
expm1-xcheck: $(EXPM1_GEN)
	common/tools/xcheck.sh expm1 $(PREFIX) $(PREFIX)
expm1-hashes: $(EXPM1_B)/hwy_test
	$(RUN) ./$(EXPM1_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += expm1-all
TEST_TARGETS   += expm1-test
HWY_TARGETS    += expm1-hwy
HASHES_TARGETS += expm1-hashes
XCHECK_TARGETS += expm1-xcheck

ifeq ($(ARCH),x86_64)
expm1-insns: $(EXPM1_B)/svml_expm1_d_ha.o $(EXPM1_B)/svml_asm.o
	common/tools/compare_insns.sh $(EXPM1_B)/svml_expm1_d_ha.o svml_expm18_ha $(EXPM1_B)/svml_asm.o __svml_expm18_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway targets vs the translation
expm1-avx512-test: $(EXPM1_B)/test_expm18 $(EXPM1_B)/hwy_test
	$(RUN) ./$(EXPM1_B)/test_expm18 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(EXPM1_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += expm1-insns
AVX512_TARGETS += expm1-avx512-test
endif

.PHONY: expm1-all expm1-test expm1-hwy expm1-hashes expm1-xcheck expm1-insns expm1-avx512-test
