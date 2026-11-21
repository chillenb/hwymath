# hwymath

[![CI](https://github.com/chillenb/hwymath/actions/workflows/ci.yml/badge.svg)](https://github.com/chillenb/hwymath/actions/workflows/ci.yml)

Studies of Intel's open-source high-accuracy SVML math functions
([numpy/SVML](https://github.com/numpy/SVML), AVX-512 assembly): for each function, a bit-identical
translation to C, ports to other instruction sets, improved variants, and portable
[Google Highway](https://github.com/google/highway) versions, each tested bit for bit against Intel's
code and compared with Intel's libsvml, glibc, SLEEF and Highway's own math library.

| functions (double) | versions | highlights |
|---|---|---|
| [exp](functions/exp/README.md) | C translation, AVX2 ports and variants, Highway hp | more accurate than Intel's (max 0.505–0.511 ulp vs 0.80), same speed on AVX2 |
| [exp2](functions/exp2/README.md), [expm1](functions/expm1/README.md) | C translation, Highway | bit-identical to Intel's AVX-512 code on every FMA target |
| [log](functions/log/README.md), [log2](functions/log2/README.md), [log10](functions/log10/README.md), [log1p](functions/log1p/README.md) | C translation, Highway | exact model of the rounded `vrcp14pd`, verified against Intel's reference algorithm |
| [pow](functions/pow/README.md) | C translation, Highway | exact `{rz-sae}`/`{rd-sae}` arithmetic and `vreducepd`; checked per operation against the instructions |
| [cbrt](functions/cbrt/README.md) | C translation, Highway | bit-identical; exact model of the rounded `vrcp14pd`, verified against the instruction |
| [sin](functions/sin/README.md), [cos](functions/cos/README.md), [tan](functions/tan/README.md) | C translation, Highway | Payne–Hanek reduction for huge arguments; exact `vrcp14pd` for tan |
| [asin](functions/asin/README.md), [acos](functions/acos/README.md), [atan](functions/atan/README.md), [atan2](functions/atan2/README.md) | C translation (callouts translated mechanically), Highway | exact `vrsqrt14pd`/`vrcp14pd` emulation after Intel's `recip14.c` |
| [sinh](functions/sinh/README.md), [cosh](functions/cosh/README.md), [tanh](functions/tanh/README.md) | C translation, Highway | as fast as or faster than Intel's HA code on AVX-512; sinh and cosh also on AVX2 |
| [asinh](functions/asinh/README.md), [acosh](functions/acosh/README.md), [atanh](functions/atanh/README.md) | C translation, Highway | bit-identical, as fast as Intel's on AVX-512 |

Every C translation is bit-identical to Intel's assembly (vector code and scalar callout), and every Highway
version of an SVML algorithm is bit-identical to Intel's AVX-512 code on every target with FMA (AVX-512,
AVX2, NEON, SVE and SVE2 at any vector length), checked in CI on x86-64 and AArch64 hardware.

All 67 high-accuracy sources and their status: [docs/svml-ha-sources.md](docs/svml-ha-sources.md).
Accuracy and speed: [Performance.md](Performance.md).

## Getting started

Everything (compilers, Highway, SLEEF, Intel's runtime for comparisons) comes from conda-forge via
[pixi](https://pixi.sh); Intel's assembly is the `third_party/SVML` submodule.

```bash
git clone --recursive https://github.com/chillenb/hwymath && cd hwymath
pixi run make test        # correctness checks of every function
pixi run make hwy         # Highway versions on every supported target, with comparisons and benchmarks
pixi run make compare     # x86-64: accuracy comparisons with Intel SVML, glibc, libmvec, SLEEF
pixi run make bench       # x86-64: throughput
pixi run make exp-test    # one function only: <function>-<target>
pixi run -e cross make xcheck   # AArch64 NEON/SVE/SVE2 under QEMU vs x86-64, bit for bit
```

See the comment at the top of the `Makefile` for all targets.

## Layout

```
functions/<name>/        one directory per SVML function (all precisions studied for it)
  README.md              algorithm, versions, tests, notes
  Performance.md         accuracy and speed of all versions and of the libraries compared with
  function.mk            build rules and <name>-* make targets (included by the top-level Makefile)
  svml/                  bit-identical C translations of svml_z0_<name>_<d|s|h>_ha.s + constants
  avx2/, hwy/, ...       ports and variants
  tests/, tools/         tests, benchmarks, comparison implementations; generators
common/                  code used across functions (see common/README.md)
  include/hwymath/       test helpers (bit casts, RNG, accuracy statistics, hashes, timing); hwy/: exact
                         Highway emulations of AVX-512 instructions (vrcp14pd, vrsqrt14pd, ...)
  mk/                    shared build settings (config.mk) and rules for the shared code (common.mk)
  tests/                 tests of the shared code
  tools/                 xcheck.sh (AArch64 under QEMU vs x86-64), extract_tables.py, compare_insns.sh,
                         check_hashes.py (CI cross-machine check), svml_catalog.py
docs/porting.md          how a function is ported and verified (also the brief for parallel agents)
docs/svml-ha-sources.md  catalog of all SVML high-accuracy sources (make catalog)
third_party/SVML/        numpy/SVML submodule
```

### Adding a function

1. Create `functions/<name>/` with a `function.mk` modelled on `functions/exp/function.mk`: build into
   `build/<name>/`, define `<name>-all`, `<name>-test`, ... and append them to the aggregate variables
   (`ALL_TARGETS`, `TEST_TARGETS`, `HWY_TARGETS`, `HASHES_TARGETS`, and on x86-64 `COMPARE_TARGETS`,
   `BENCH_TARGETS`, `INSNS_TARGETS`, `AVX512_TARGETS`, ...). The top-level Makefile picks it up
   automatically.
2. Translate `svml_z0_<name>_<p>_ha.s` to `functions/<name>/svml/svml_<name>_<p>_ha.c`; extract its
   constants with `common/tools/extract_tables.py`; check the instruction mix with
   `common/tools/compare_insns.sh`.
3. Make `<name>-hashes` append `arch <name> target variant n hash` lines to `$(HASH_FILE)`; CI then
   checks that results agree across machines.
4. `make catalog` to update `docs/svml-ha-sources.md`; add a row to the table above.

## CI

`.github/workflows/ci.yml` runs `make test`, `make hashes`, `make hwy`, `make compare`, `make bench` and
`make insns` on GitHub's x86-64 and AArch64 runners, for pull requests and pushes to `main`. AVX-512 runs natively when the runner has it
(usually); otherwise `make avx512-test` runs under Intel SDE. A final job checks that the result hashes
of every function agree across machines and targets.

## License

BSD-3-Clause (see `LICENSE`). The algorithms, constants and scalar routines are derived from Intel's
SVML (Copyright (C) 2023 Intel Corporation, BSD-3-Clause); the `vrcp14pd`/`vrsqrt14pd` coefficient tables come
from Intel's `recip14.c` (Copyright (c) 2015 Intel Corporation, BSD-3-Clause).
