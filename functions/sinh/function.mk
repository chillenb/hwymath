# sinh, double precision: bit-identical C/AVX-512 translation of svml_z0_sinh_d_ha.s and a Google
# Highway version of the same algorithm. Included by the top-level Makefile; paths are relative to the
# repository root. Targets: sinh-all sinh-test sinh-hwy sinh-hashes sinh-xcheck (AArch64 under QEMU, in the
# pixi "cross" environment), and on x86-64 also sinh-insns sinh-avx512-test.

SINH     := functions/sinh
SINH_B   := build/sinh
SINH_ASM := $(SVML_ASM)/svml_z0_sinh_d_ha.s
SINH_H   := $(SINH)/sinh_d.h
SINH_GEN := $(SINH)/svml/svml_sinh_d_ha_constants.inc

# generated constants (committed; regenerated when the generator changes)
$(SINH)/svml/svml_sinh_d_ha_constants.inc: $(SINH)/tools/gen_svml_constants.py common/tools/extract_tables.py $(SINH_ASM)
	python3 -I $< $(SINH_ASM) > $@

ifeq ($(ARCH),x86_64)
SINH_OBJS  := $(SINH_B)/svml_sinh_d_ha.o $(SINH_B)/svml_asm.o
SINH_TESTS := $(SINH_B)/test_rare $(SINH_B)/test_sinh8
else
SINH_OBJS  := $(SINH_B)/svml_sinh_d_ha.o
SINH_TESTS :=
endif
SINH_PROGS := $(SINH_TESTS) $(SINH_B)/hwy_test

$(SINH_B):
	mkdir -p $@

$(SINH_B)/svml_sinh_d_ha.o: $(SINH)/svml/svml_sinh_d_ha.c $(SINH_H) $(SINH_GEN) | $(SINH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(SINH_B)/svml_asm.o: $(SINH_ASM) | $(SINH_B)
	sed 's/^__svml_dsinh_ha_cout_rare_internal:/\t.globl __svml_dsinh_ha_cout_rare_internal\n&/' $< > $(SINH_B)/svml_asm.s
	$(CC) -c $(SINH_B)/svml_asm.s -o $@

SINH_TEST_DEPS := $(SINH)/tests/inputs.h common/include/hwymath/test_util.h $(SINH_OBJS)
$(SINH_B)/test_%: $(SINH)/tests/test_%.c $(SINH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(SINH_OBJS) -lm -o $@

# Highway version
$(SINH_B)/hwy_sinh.o: $(SINH)/hwy/hwy_sinh.cc $(SINH)/hwy/hwy_sinh_constants.h $(SINH_GEN) | $(SINH_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(SINH_B)/hwy_test: $(SINH)/tests/hwy_test.cc $(SINH_TEST_DEPS) common/include/hwymath/accuracy.h $(SINH_B)/hwy_sinh.o
	$(CXX) $(CXXFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) $< $(SINH_B)/hwy_sinh.o $(SINH_OBJS) $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

sinh-all: $(SINH_PROGS)
sinh-test: $(SINH_TESTS) $(SINH_B)/hwy_test
	@set -e; for t in $(SINH_B)/test_rare; do [ -e $$t ] || continue; $(RUN) ./$$t; done
	$(RUN) ./$(SINH_B)/hwy_test 2000000 --no-bench
sinh-hwy: $(SINH_B)/hwy_test
	./$(SINH_B)/hwy_test
sinh-hashes: $(SINH_B)/hwy_test
	$(RUN) ./$(SINH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
sinh-xcheck: $(SINH_GEN)
	common/tools/xcheck.sh sinh $(PREFIX) $(PREFIX)
ALL_TARGETS    += sinh-all
TEST_TARGETS   += sinh-test
HWY_TARGETS    += sinh-hwy
HASHES_TARGETS += sinh-hashes
XCHECK_TARGETS += sinh-xcheck

ifeq ($(ARCH),x86_64)
sinh-insns: $(SINH_B)/svml_sinh_d_ha.o $(SINH_B)/svml_asm.o
	common/tools/compare_insns.sh $(SINH_B)/svml_sinh_d_ha.o svml_sinh8_ha $(SINH_B)/svml_asm.o __svml_sinh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets vs the translation
sinh-avx512-test: $(SINH_B)/test_sinh8 $(SINH_B)/hwy_test
	$(RUN) ./$(SINH_B)/test_sinh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(SINH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += sinh-insns
AVX512_TARGETS += sinh-avx512-test
endif

.PHONY: sinh-all sinh-test sinh-hwy sinh-hashes sinh-xcheck sinh-insns sinh-avx512-test
