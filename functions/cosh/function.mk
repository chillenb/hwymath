# cosh, double precision: bit-identical C/AVX-512 translation of svml_z0_cosh_d_ha.s and a Google
# Highway version of the same algorithm. Included by the top-level Makefile; paths are relative to the
# repository root. Targets: cosh-all cosh-test cosh-hwy cosh-hashes cosh-xcheck (AArch64 under QEMU, in the
# pixi "cross" environment), and on x86-64 also cosh-insns cosh-avx512-test.

COSH     := functions/cosh
COSH_B   := build/cosh
COSH_ASM := $(SVML_ASM)/svml_z0_cosh_d_ha.s
COSH_H   := $(COSH)/cosh_d.h
COSH_GEN := $(COSH)/svml/svml_cosh_d_ha_constants.inc

# generated constants (committed; regenerated when the generator changes)
$(COSH)/svml/svml_cosh_d_ha_constants.inc: $(COSH)/tools/gen_svml_constants.py common/tools/extract_tables.py $(COSH_ASM)
	python3 -I $< $(COSH_ASM) > $@

ifeq ($(ARCH),x86_64)
COSH_OBJS  := $(COSH_B)/svml_cosh_d_ha.o $(COSH_B)/svml_asm.o
COSH_TESTS := $(COSH_B)/test_rare $(COSH_B)/test_cosh8
else
COSH_OBJS  := $(COSH_B)/svml_cosh_d_ha.o
COSH_TESTS :=
endif
COSH_PROGS := $(COSH_TESTS) $(COSH_B)/hwy_test

$(COSH_B):
	mkdir -p $@

$(COSH_B)/svml_cosh_d_ha.o: $(COSH)/svml/svml_cosh_d_ha.c $(COSH_H) $(COSH_GEN) | $(COSH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(COSH_B)/svml_asm.o: $(COSH_ASM) | $(COSH_B)
	sed 's/^__svml_dcosh_ha_cout_rare_internal:/\t.globl __svml_dcosh_ha_cout_rare_internal\n&/' $< > $(COSH_B)/svml_asm.s
	$(CC) -c $(COSH_B)/svml_asm.s -o $@

COSH_TEST_DEPS := $(COSH)/tests/inputs.h common/include/hwymath/test_util.h $(COSH_OBJS)
$(COSH_B)/test_%: $(COSH)/tests/test_%.c $(COSH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(COSH_OBJS) -lm -o $@

# Highway version
$(COSH_B)/hwy_cosh.o: $(COSH)/hwy/hwy_cosh.cc $(COSH)/hwy/hwy_cosh_constants.h $(COSH_GEN) | $(COSH_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(COSH_B)/hwy_test: $(COSH)/tests/hwy_test.cc $(COSH_TEST_DEPS) common/include/hwymath/accuracy.h $(COSH_B)/hwy_cosh.o
	$(CXX) $(CXXFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) $< $(COSH_B)/hwy_cosh.o $(COSH_OBJS) $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

cosh-all: $(COSH_PROGS)
cosh-test: $(COSH_TESTS) $(COSH_B)/hwy_test
	@set -e; for t in $(COSH_B)/test_rare; do [ -e $$t ] || continue; $(RUN) ./$$t; done
	$(RUN) ./$(COSH_B)/hwy_test 2000000 --no-bench
cosh-hwy: $(COSH_B)/hwy_test
	./$(COSH_B)/hwy_test
cosh-hashes: $(COSH_B)/hwy_test
	$(RUN) ./$(COSH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
cosh-xcheck: $(COSH_GEN)
	./$(COSH)/tools/xcheck.sh $(PREFIX) $(PREFIX)
ALL_TARGETS    += cosh-all
TEST_TARGETS   += cosh-test
HWY_TARGETS    += cosh-hwy
HASHES_TARGETS += cosh-hashes
XCHECK_TARGETS += cosh-xcheck

ifeq ($(ARCH),x86_64)
cosh-insns: $(COSH_B)/svml_cosh_d_ha.o $(COSH_B)/svml_asm.o
	common/tools/compare_insns.sh $(COSH_B)/svml_cosh_d_ha.o svml_cosh8_ha $(COSH_B)/svml_asm.o __svml_cosh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets vs the translation
cosh-avx512-test: $(COSH_B)/test_cosh8 $(COSH_B)/hwy_test
	$(RUN) ./$(COSH_B)/test_cosh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(COSH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += cosh-insns
AVX512_TARGETS += cosh-avx512-test
endif

.PHONY: cosh-all cosh-test cosh-hwy cosh-hashes cosh-xcheck cosh-insns cosh-avx512-test
