# atan, double precision: bit-identical C/AVX-512 translation of svml_z0_atan_d_ha.s and a Google Highway
# version. Included by the top-level Makefile; paths are relative to the repository root.
# Targets: atan-all atan-test atan-hwy atan-hashes, and on x86-64 also atan-insns atan-avx512-test.
# The scalar callout is translated by functions/asin/tools/asm2c.py (machine model functions/asin/svml/asm_model.h).

ATAN     := functions/atan
ATAN_B   := build/atan
ATAN_ASM := $(SVML_ASM)/svml_z0_atan_d_ha.s
ATAN_H   := $(ATAN)/atan_d.h
ATAN_ASM2C  := functions/asin/tools/asm2c.py
ATAN_MODEL_H := functions/asin/svml/asm_model.h

# --- generated sources (committed; regenerated when their generators change) ---
ATAN_GEN := $(ATAN)/svml/svml_atan_d_ha_constants.inc $(ATAN)/svml/svml_datan_ha_cout_rare.inc
$(ATAN)/svml/svml_atan_d_ha_constants.inc: $(ATAN)/tools/gen_svml_constants.py common/tools/extract_tables.py $(ATAN_ASM)
	python3 -I $< $(ATAN_ASM) > $@
$(ATAN)/svml/svml_datan_ha_cout_rare.inc: $(ATAN_ASM2C) common/tools/extract_tables.py $(ATAN_ASM)
	python3 -I $(ATAN_ASM2C) $(ATAN_ASM) __svml_datan_ha_cout_rare_internal svml_datan_ha_cout_rare 2 > $@

# --- what gets built on this architecture ---
ATAN_OBJS := $(ATAN_B)/svml_atan_d_ha.o
ifeq ($(ARCH),x86_64)
ATAN_OBJS  += $(ATAN_B)/svml_asm.o
ATAN_TESTS := $(ATAN_B)/test_rare $(ATAN_B)/test_atan8
else
ATAN_TESTS :=
endif
ATAN_PROGS := $(ATAN_TESTS) $(ATAN_B)/hwy_test

$(ATAN_B):
	mkdir -p $@

$(ATAN_B)/svml_atan_d_ha.o: $(ATAN)/svml/svml_atan_d_ha.c $(ATAN_H) $(ATAN_MODEL_H) $(ATAN_GEN) | $(ATAN_B)
	$(CC) $(CFLAGS) -I$(dir $(ATAN_MODEL_H)) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
$(ATAN_B)/svml_asm.o: $(ATAN_ASM) | $(ATAN_B)
	sed 's/^__svml_datan_ha_cout_rare_internal:/\t.globl __svml_datan_ha_cout_rare_internal\n&/' $< > $(ATAN_B)/svml_asm.s
	$(CC) -c $(ATAN_B)/svml_asm.s -o $@

ATAN_TEST_DEPS := $(ATAN)/tests/inputs.h common/include/hwymath/test_util.h $(ATAN_OBJS)
$(ATAN_B)/test_%: $(ATAN)/tests/test_%.c $(ATAN_TEST_DEPS)
	$(CC) $(CFLAGS) $(INC) $< $(ATAN_OBJS) -lm -o $@

# --- Highway version ---
$(ATAN_B)/hwy_atan.o: $(ATAN)/hwy/hwy_atan.cc $(ATAN)/hwy/hwy_atan_constants.h $(ATAN_GEN) | $(ATAN_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(ATAN_B)/extras.o: $(ATAN)/tests/extras.c $(ATAN)/tests/extras.h | $(ATAN_B)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
$(ATAN_B)/hwy_test: $(ATAN)/tests/hwy_test.cc $(ATAN)/tests/extras.h common/include/hwymath/accuracy.h \
                    $(ATAN_TEST_DEPS) $(ATAN_B)/hwy_atan.o $(ATAN_B)/extras.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(ATAN_B)/hwy_atan.o $(ATAN_OBJS) $(ATAN_B)/extras.o $(LIBPATH) -lhwy \
	  $(SLEEF_LIB) $(INTEL_LIB) -lm -o $@

# --- entry points ---
atan-all: $(ATAN_PROGS)
atan-test: $(ATAN_TESTS) $(ATAN_B)/hwy_test
	@set -e; for t in $(ATAN_TESTS); do $(RUN) ./$$t || [ $$? -eq 77 ]; done
	$(RUN) ./$(ATAN_B)/hwy_test 2000000 --no-bench
atan-hwy: $(ATAN_B)/hwy_test
	./$(ATAN_B)/hwy_test
atan-hashes: $(ATAN_B)/hwy_test
	$(RUN) ./$(ATAN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
ALL_TARGETS    += atan-all
TEST_TARGETS   += atan-test
HWY_TARGETS    += atan-hwy
HASHES_TARGETS += atan-hashes

ifeq ($(ARCH),x86_64)
atan-insns: $(ATAN_B)/svml_atan_d_ha.o $(ATAN_B)/svml_asm.o
	common/tools/compare_insns.sh $(ATAN_B)/svml_atan_d_ha.o svml_atan8_ha $(ATAN_B)/svml_asm.o __svml_atan8_ha
# AVX-512: bit-exactness of the translation vs Intel's assembly, Highway AVX-512 targets
atan-avx512-test: $(ATAN_B)/test_atan8 $(ATAN_B)/hwy_test
	$(RUN) ./$(ATAN_B)/test_atan8 2000000 || [ $$? -eq 77 ]
	$(RUN) ./$(ATAN_B)/hwy_test $(HASH_INPUTS) --no-bench --hash $(HASH_FILE)
INSNS_TARGETS  += atan-insns
AVX512_TARGETS += atan-avx512-test
endif

.PHONY: atan-all atan-test atan-hwy atan-hashes atan-insns atan-avx512-test
