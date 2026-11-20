# Settings shared by all functions; included by the top-level Makefile. Run make inside pixi
# (`pixi run make ...`): compilers, Highway, SLEEF and Intel's libsvml come from $(PREFIX).

ARCH     := $(shell uname -m)
PREFIX   ?= $(CONDA_PREFIX)
SVML_DIR ?= third_party/SVML
SVML_ASM := $(SVML_DIR)/linux/avx512

# -ffp-contract=off: scalar reference code relies on unfused operations, and test inputs must be
# bit-identical on every architecture
CFLAGS   := -O2 -g -Wall -Wextra -ffp-contract=off
CXXFLAGS := -O2 -g -std=c++17 -Wall -Wextra -ffp-contract=off
INC      := -I. -Icommon/include -I$(PREFIX)/include
LIBPATH  := -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib

# optional comparison libraries
ifneq ($(wildcard $(PREFIX)/include/sleef.h),)
SLEEF_DEF := -DHAVE_SLEEF
SLEEF_LIB := -lsleef
endif
ifneq ($(wildcard $(PREFIX)/lib/libsvml.so),)
INTEL_DEF := -DHAVE_INTEL_SVML
INTEL_LIB := -lsvml -lintlc
endif

# prefix for running test programs, e.g. RUN="$SDE_PATH/sde64 -spr --" for AVX-512 under Intel SDE
RUN ?=
# bit-identity hashes (make hashes): number of random inputs and output file
HASH_INPUTS ?= 2000000
HASH_FILE   ?= hashes-$(ARCH).txt
