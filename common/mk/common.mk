# Code shared by several functions (common/include/hwymath, common/tests). Included by the top-level Makefile
# before the functions' function.mk files. Targets: common-all common-test, and on x86-64 also
# common-avx512-test.

COMMON_B := build/common

# exact Highway emulations of AVX-512 instructions (VRCP14PD, VRSQRT14PD, VGETEXPPD, VPERMT2PD) with their
# generated data; on AVX-512 targets the header uses the instructions themselves
AVX512_EMUL := common/include/hwymath/hwy/avx512_emul-inl.h common/include/hwymath/approx14_tables.h

common/include/hwymath/approx14_tables.h: common/tools/gen_approx14.py
	python3 -I $< > $@

$(COMMON_B):
	mkdir -p $@

# test_approx14: Intel's reference algorithm vs the instructions (with AVX-512 or under SDE), and the Highway
# emulation vs the reference on every target, for all 2^16 mantissa prefixes
$(COMMON_B)/hwy_approx14.o: common/tests/approx14/hwy_approx14.cc $(AVX512_EMUL) | $(COMMON_B)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
$(COMMON_B)/test_approx14: common/tests/approx14/test_approx14.cc common/tests/approx14/approx14_ref.h \
                           $(COMMON_B)/hwy_approx14.o
	$(CXX) $(CXXFLAGS) $(INC) $< $(COMMON_B)/hwy_approx14.o $(LIBPATH) -lhwy -lm -o $@

common-all: $(COMMON_B)/test_approx14
common-test: $(COMMON_B)/test_approx14
	$(RUN) ./$(COMMON_B)/test_approx14
ALL_TARGETS  += common-all
TEST_TARGETS += common-test

ifeq ($(ARCH),x86_64)
common-avx512-test: common-test
AVX512_TARGETS += common-avx512-test
endif

.PHONY: common-all common-test common-avx512-test
