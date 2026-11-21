# tanh, double precision: bit-identical C/AVX-512 translation of svml_z0_tanh_d_ha.s and a Google
# Highway version of the same algorithm. Included by the top-level Makefile; paths are relative to the
# repository root. Targets: tanh-all tanh-test tanh-hwy tanh-hashes tanh-xcheck (AArch64 under QEMU, in the
# pixi "cross" environment), and on x86-64 also tanh-insns tanh-avx512-test.

TANH     := functions/tanh
TANH_B   := build/tanh
TANH_ASM := $(SVML_ASM)/svml_z0_tanh_d_ha.s
TANH_H   := $(TANH)/tanh_d.h
TANH_GEN := $(TANH)/svml/svml_tanh_d_ha_constants.inc

# generated constants (committed; regenerated when the generator changes)
$(TANH)/svml/svml_tanh_d_ha_constants.inc: $(TANH)/tools/gen_svml_constants.py common/tools/extract_tables.py $(TANH_ASM)
	python3 -I $< $(TANH_ASM) > $@

ifeq ($(ARCH),x86_64)
TANH_OBJS  := $(TANH_B)/svml_tanh_d_ha.o $(TANH_B)/svml_asm.o
TANH_TESTS := $(TANH_B)/test_rare $(TANH_B)/test_tanh8
else
TANH_OBJS  := $(TANH_B)/svml_tanh_d_ha.o
TANH_TESTS :=
endif
TANH_PROGS := $(TANH_TESTS) $(TANH_B)/hwy_test

$(TANH_B):
	mkdir -p $@

$(TANH_B)/svml_tanh_d_ha.o: $(TANH)/svml/svml_tanh_d_ha.c $(TANH_H) $(TANH_GEN) | $(TANH_B)
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(TANH_B)/svml_asm.o: $(TANH_ASM) | $(TANH_B)
	sed 's/^__svml_dtanh_ha_cout_rare_internal:/\t.globl __svml_dtanh_ha_cout_rare_internal\n&/' $< > $(TANH_B)/svml_asm.s
	$(CC) -c $(TANH_B)/svml_asm.s -o $@

TANH_TEST_DEPS := $(TANH)/tests/inputs.h common/include/hwymath/test_util.h $(TANH_OBJS)
$(TANH_B)/test_%: $(TANH)/tests/test_%.c $(TANH_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(TANH_OBJS) -lm -o $@

# Highway version
$(TANH_B)/hwy_tanh.o: $(TANH)/hwy/hwy_tanh.cc $(TANH)/hwy/hwy_tanh_constants.h $(TANH_GEN) | $(TANH_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(TANH_B)/hwy_test: $(TANH)/tests/hwy_test.cc $(TANH_TEST_DEPS) common/include/hwymath/accuracy.h $(TANH_B)/hwy_tanh.o
	$(CXX) $(CXXFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) $< $(TANH_B)/hwy_tanh.o $(TANH_OBJS) $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

tanh-all: $(TANH_PROGS)
tanh-test: $(TANH_TESTS) $(TANH_B)/hwy_test
	@set -e; for t in $(TANH_B)/test_rare; do [ -e $$t ] || continue; $(RUN) ./$$t; done
	$(RUN) ./$(TANH_B)/hwy_test 2000000 --no-bench
tanh-hwy: $(TANH_B)/hwy_test
	./$(TANH_B)/hwy_test
tanh-hashes: $(TANH_B)/hwy_test
	$(RUN) ./$(TANH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
tanh-xcheck: $(TANH_GEN)
	common/tools/xcheck.sh tanh $(PREFIX) $(PREFIX)
ALL_TARGETS    += tanh-all
TEST_TARGETS   += tanh-test
HWY_TARGETS    += tanh-hwy
HASHES_TARGETS += tanh-hashes
XCHECK_TARGETS += tanh-xcheck

ifeq ($(ARCH),x86_64)
tanh-insns: $(TANH_B)/svml_tanh_d_ha.o $(TANH_B)/svml_asm.o
	common/tools/compare_insns.sh $(TANH_B)/svml_tanh_d_ha.o svml_tanh8_ha $(TANH_B)/svml_asm.o __svml_tanh8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets vs the translation
tanh-avx512-test: $(TANH_B)/test_tanh8 $(TANH_B)/hwy_test
	$(RUN) ./$(TANH_B)/test_tanh8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(TANH_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += tanh-insns
AVX512_TARGETS += tanh-avx512-test
endif

.PHONY: tanh-all tanh-test tanh-hwy tanh-hashes tanh-xcheck tanh-insns tanh-avx512-test
