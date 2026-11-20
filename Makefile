# hwymath: studies of Intel SVML's high-accuracy math functions (one directory per function under
# functions/), with C translations, AVX2 ports and portable Google Highway versions.
#
# Run inside the pixi environment:
#   pixi run make test         correctness checks of every function (exit status != 0 on failure)
#   pixi run make hwy          Highway versions on every compiled+supported target, with comparisons
#   pixi run make compare      x86-64: accuracy comparisons with Intel SVML, glibc, libmvec, SLEEF
#   pixi run make bench        x86-64: throughput
#   pixi run make insns        x86-64: instruction mix of the C translations vs Intel's assembly
#   pixi run make hashes       append result hashes to $(HASH_FILE) (CI cross-machine check)
#   pixi run make avx512-test  AVX-512 tests only (CI runs them under Intel SDE without AVX-512)
#   pixi run -e cross make xcheck   Highway versions on AArch64 NEON/SVE/SVE2 under QEMU vs x86-64
#   pixi run make catalog      regenerate docs/svml-ha-sources.md
# Per function: <function>-<target>, e.g. make exp-test.

all:

include common/mk/config.mk

FUNCTIONS := $(sort $(notdir $(patsubst %/function.mk,%,$(wildcard functions/*/function.mk))))

# aggregate targets; each functions/<f>/function.mk appends its own <f>-* targets
ALL_TARGETS :=
TEST_TARGETS :=
HWY_TARGETS :=
COMPARE_TARGETS :=
BENCH_TARGETS :=
INSNS_TARGETS :=
HASHES_TARGETS :=
AVX512_TARGETS :=
XCHECK_TARGETS :=
include $(foreach f,$(FUNCTIONS),functions/$(f)/function.mk)

all: $(ALL_TARGETS)
test: $(TEST_TARGETS)
hwy: $(HWY_TARGETS)
compare: $(COMPARE_TARGETS)
bench: $(BENCH_TARGETS)
insns: $(INSNS_TARGETS)
hashes: $(HASHES_TARGETS)
avx512-test: $(AVX512_TARGETS)
xcheck: $(XCHECK_TARGETS)

catalog:
	python3 -I common/tools/svml_catalog.py > docs/svml-ha-sources.md

clean:
	rm -rf build

.PHONY: all test hwy compare bench insns hashes avx512-test xcheck catalog clean
