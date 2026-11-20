# acos, double precision: bit-identical C/AVX-512 translation of svml_z0_acos_d_ha.s and a Google Highway
# version. Included by the top-level Makefile; paths are relative to the repository root.
# Targets: acos-all acos-test acos-hwy acos-hashes, and on x86-64 also acos-insns acos-avx512-test.
# The scalar callout is translated by functions/asin/tools/asm2c.py (machine model functions/asin/svml/asm_model.h).

ACOS     := functions/acos
ACOS_B   := build/acos
ACOS_ASM := $(SVML_ASM)/svml_z0_acos_d_ha.s
ACOS_H   := $(ACOS)/acos_d.h
ACOS_ASM2C  := functions/asin/tools/asm2c.py
ACOS_MODEL_H := functions/asin/svml/asm_model.h

# --- generated sources (committed; regenerated when their generators change) ---
ACOS_GEN := $(ACOS)/svml/svml_acos_d_ha_constants.inc $(ACOS)/svml/svml_dacos_ha_cout_rare.inc
$(ACOS)/svml/svml_acos_d_ha_constants.inc: $(ACOS)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ACOS_ASM)
	python3 -I $< $(ACOS_ASM) > $@
$(ACOS)/svml/svml_dacos_ha_cout_rare.inc: $(ACOS_ASM2C) common/tools/extract_tables.py $(ACOS_ASM)
	python3 -I $(ACOS_ASM2C) $(ACOS_ASM) __svml_dacos_ha_cout_rare_internal svml_dacos_ha_cout_rare 2 > $@

# --- what gets built on this architecture ---
ACOS_OBJS := $(ACOS_B)/svml_acos_d_ha.o
ifeq ($(ARCH),x86_64)
ACOS_OBJS  += $(ACOS_B)/svml_asm.o
ACOS_TESTS := $(ACOS_B)/test_rare $(ACOS_B)/test_acos8
else
ACOS_TESTS :=
endif
ACOS_PROGS := $(ACOS_TESTS) $(ACOS_B)/hwy_test

$(ACOS_B):
	mkdir -p $@

$(ACOS_B)/svml_acos_d_ha.o: $(ACOS)/svml/svml_acos_d_ha.c $(ACOS_H) $(ACOS_MODEL_H) $(ACOS_GEN) | $(ACOS_B)
	$(CC) $(CFLAGS) -I$(dir $(ACOS_MODEL_H)) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ACOS_B)/svml_asm.o: $(ACOS_ASM) | $(ACOS_B)
	sed 's/^__svml_dacos_ha_cout_rare_internal:/\t.globl __svml_dacos_ha_cout_rare_internal\n&/' $< > $(ACOS_B)/svml_asm.s
	$(CC) -c $(ACOS_B)/svml_asm.s -o $@

ACOS_TEST_DEPS := $(ACOS)/tests/inputs.h common/include/hwymath/test_util.h $(ACOS_OBJS)
$(ACOS_B)/test_%: $(ACOS)/tests/test_%.c $(ACOS_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ACOS_OBJS) -lm -o $@

# --- Highway version ---
$(ACOS_B)/hwy_acos.o: $(ACOS)/hwy/hwy_acos.cc $(ACOS)/hwy/hwy_acos_constants.h $(ACOS_GEN) | $(ACOS_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(ACOS_B)/extras.o: $(ACOS)/tests/extras.c $(ACOS)/tests/extras.h | $(ACOS_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ACOS_B)/hwy_test: $(ACOS)/tests/hwy_test.cc $(ACOS)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(ACOS_TEST_DEPS) $(ACOS_B)/hwy_acos.o $(ACOS_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(ACOS_B)/hwy_acos.o $(ACOS_OBJS) $(ACOS_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

# --- entry points ---
acos-all: $(ACOS_PROGS)
acos-test: $(ACOS_TESTS) $(ACOS_B)/hwy_test
	@set -e; for t in $(ACOS_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ACOS_B)/hwy_test 2000000 --no-bench
acos-hwy: $(ACOS_B)/hwy_test
	./$(ACOS_B)/hwy_test
acos-hashes: $(ACOS_B)/hwy_test
	$(RUN) ./$(ACOS_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += acos-all
TEST_TARGETS   += acos-test
HWY_TARGETS    += acos-hwy
HASHES_TARGETS += acos-hashes

ifeq ($(ARCH),x86_64)
acos-insns: $(ACOS_B)/svml_acos_d_ha.o $(ACOS_B)/svml_asm.o
	common/tools/compare_insns.sh $(ACOS_B)/svml_acos_d_ha.o svml_acos8_ha $(ACOS_B)/svml_asm.o __svml_acos8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
acos-avx512-test: $(ACOS_B)/test_acos8 $(ACOS_B)/hwy_test
	$(RUN) ./$(ACOS_B)/test_acos8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ACOS_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += acos-insns
AVX512_TARGETS += acos-avx512-test
endif

.PHONY: acos-all acos-test acos-hwy acos-hashes acos-insns acos-avx512-test
