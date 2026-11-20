# hwymath: SVML double exp in C, AVX2 ports, and a portable Google Highway version.
# Run inside the pixi environment (compilers, Highway, SLEEF and Intel's libsvml come from conda-forge):
#   pixi run make test      bit-exactness and accuracy checks (exit status != 0 on failure)
#   pixi run make hwy       Highway exp on every compiled+supported target, with benchmarks and comparisons
#   pixi run make compare   x86-64: AVX2 variants vs Intel SVML, glibc, libmvec, SLEEF (accuracy)
#   pixi run make bench     x86-64: throughput of the same
#   pixi run make insns     x86-64: instruction mix of the AVX-512 translation vs Intel's assembly
#   pixi run -e cross make xcheck   Highway exp on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64, bit for bit

ARCH     := $(shell uname -m)
PREFIX   ?= $(CONDA_PREFIX)
SVML_DIR ?= third_party/SVML
ASM      := $(SVML_DIR)/linux/avx512/svml_z0_exp_d_ha.s
# -ffp-contract=off: the scalar code relies on unfused operations, and test inputs must be identical everywhere
CFLAGS   := -O2 -g -Wall -Wextra -ffp-contract=off
CXXFLAGS := -O2 -g -std=c++17 -Wall -Wextra -ffp-contract=off
INC      := -I. -I$(PREFIX)/include
LIBPATH  := -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib

ifneq ($(wildcard $(PREFIX)/include/sleef.h),)
SLEEF_DEF := -DHAVE_SLEEF
SLEEF_LIB := -lsleef
endif
ifneq ($(wildcard $(PREFIX)/lib/libsvml.so),)
INTEL_DEF := -DHAVE_INTEL_SVML
INTEL_LIB := -lsvml -lintlc
endif

ifeq ($(ARCH),x86_64)
OBJS   := build/svml_exp_d_ha.o build/svml_exp_d_ha_avx2.o build/svml_exp_d_avx2_fast.o build/svml_exp_d_avx2_hp.o \
          build/svml_asm.o
EXTRAS := build/extras_x86.o
TESTS  := build/test_rare build/test_model build/test_exp8
PROGS  := $(TESTS) build/hwy_test $(if $(INTEL_LIB),build/compare_avx2 build/bench)
else
OBJS   := build/svml_exp_d_ha.o
EXTRAS := build/extras_arm.o build/extras_arm_sve.o
TESTS  :=
PROGS  := build/hwy_test
endif

all: $(PROGS)

build:
	mkdir -p build

# --- generated constants (committed; regenerated when their sources change) ---
constants.inc: tools/gen_constants.py tools/extract_tables.py $(ASM)
	python3 -I tools/gen_constants.py $(ASM) > $@

# fast variant: table size 2^FAST_K and polynomial degree FAST_DEG (K=2: 8, K=3: 7)
FAST_K   ?= 3
FAST_DEG ?= 7
constants_fast.inc: tools/fit_exp_k3.py
	python3 -I tools/fit_exp_k3.py $(FAST_K) $(FAST_DEG) > $@

# extra-precision variant: table size 2^HP_K, Q up to r^HP_DEG (K=2: 9, K=3: 8)
HP_K     ?= 2
HP_DEG   ?= 9
constants_hp.inc: tools/fit_exp_k3.py
	python3 -I tools/fit_exp_k3.py $(HP_K) $(HP_DEG) hp > $@

hwyexp/constants_hwy.inc: tools/fit_exp_k3.py
	( python3 -I tools/fit_exp_k3.py 2 9 hp | sed -E 's/\bh(Th|Tlr|C[0-9]+|Shifter|ShifterBase|K|Ln2Hi|Ln2Lo)\b/h2\1/g'; \
	  python3 -I tools/fit_exp_k3.py 3 8 hp | sed -E 's/\bh(Th|Tlr|C[0-9]+|Shifter|ShifterBase|K|Ln2Hi|Ln2Lo)\b/h3\1/g' ) \
	  | sed 's/^static const double \(h[23]T\)/static constexpr double \1/' > $@

# --- C kernels ---
build/%.o: %.c svml_exp_d_ha.h constants.inc constants_fast.inc constants_hp.inc | build
	$(CC) $(CFLAGS) -c $< -o $@

# Intel's assembly, with the scalar callout made global so that it can be tested directly
build/svml_asm.o: $(ASM) | build
	sed 's/^__svml_dexp_ha_cout_rare_internal:/\t.globl __svml_dexp_ha_cout_rare_internal\n&/' $< > build/svml_asm.s
	$(CC) -c build/svml_asm.s -o $@

# --- comparison implementations ---
build/extras_x86.o: tests/extras_x86.c tests/extras.h svml_exp_d_ha.h | build
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $(INTEL_DEF) -c $< -o $@
build/extras_arm.o: tests/extras_arm.c tests/extras.h | build
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -c $< -o $@
build/extras_arm_sve.o: tests/extras_arm_sve.c | build
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) -march=armv8-a+sve -c $< -o $@

# --- tests and benchmarks of the C kernels (x86-64) ---
build/test_model: tests/test_model.c tests/common.h $(OBJS)
	$(CC) $(CFLAGS) -frounding-math $< $(OBJS) -lm -o $@

build/compare_avx2 build/bench: build/%: tests/%.c tests/common.h $(OBJS)
	$(CC) $(CFLAGS) $(INC) $(SLEEF_DEF) $< $(OBJS) $(LIBPATH) $(SLEEF_LIB) $(INTEL_LIB) -lmvec -lm -o $@

build/test_%: tests/test_%.c tests/common.h $(OBJS)
	$(CC) $(CFLAGS) $< $(OBJS) -lm -o $@

# --- Highway port ---
build/hwy_exp.o: hwyexp/hwy_exp.cc hwyexp/hwy_exp_constants.h hwyexp/constants_hwy.inc constants.inc | build
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

build/hwy_test: hwyexp/hwy_test.cc tests/common.h tests/extras.h build/hwy_exp.o $(OBJS) $(EXTRAS)
	$(CXX) $(CXXFLAGS) $(INC) $< build/hwy_exp.o $(OBJS) $(EXTRAS) $(LIBPATH) -lhwy $(SLEEF_LIB) $(INTEL_LIB) \
	  $(if $(filter x86_64,$(ARCH)),-lmvec) -lm -o $@

# --- entry points ---
test: $(TESTS) build/hwy_test
	@set -e; for t in $(TESTS); do ./$$t || [ $$? -eq 77 ]; done
	./build/hwy_test 2000000 --no-bench

hwy: build/hwy_test
	./build/hwy_test

compare: build/compare_avx2
	./build/compare_avx2

bench: build/bench
	./build/bench

insns: build/svml_exp_d_ha.o build/svml_asm.o
	@./tools/compare_insns.sh build/svml_exp_d_ha.o build/svml_asm.o

xcheck: hwyexp/constants_hwy.inc | build
	./tools/xcheck.sh $(PREFIX) $(PREFIX)

clean:
	rm -rf build

.PHONY: all test hwy compare bench insns xcheck clean
