# hwymath

[![CI](https://github.com/chillenb/hwymath/actions/workflows/ci.yml/badge.svg)](https://github.com/chillenb/hwymath/actions/workflows/ci.yml)

Studies of Intel's open-source high-accuracy SVML math functions
([numpy/SVML](https://github.com/numpy/SVML), AVX-512 assembly): for each function, a bit-identical
translation to C, ports to other instruction sets, improved variants, and portable
[Google Highway](https://github.com/google/highway) versions, each tested bit for bit against Intel's
code and compared with Intel's libsvml, glibc, SLEEF and Highway's own math library.

| function | status | highlights |
|---|---|---|
| [exp](functions/exp/README.md) (double) | C translation, AVX2 ports and variants, Highway | Highway version: max 0.505–0.511 ulp, bit-identical on AVX2/NEON/SVE/SVE2, as fast as Intel's AVX2 HA exp |

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
common/
  include/hwymath/       helpers shared by all tests (bit casts, deterministic RNG)
  mk/config.mk           shared build settings (pixi prefix, flags, optional libraries, RUN prefix)
  tools/                 extract_tables.py (data tables of any SVML .s file), compare_insns.sh,
                         check_hashes.py (CI cross-machine check), svml_catalog.py
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
`make insns` on GitHub's x86-64 and AArch64 runners. AVX-512 runs natively when the runner has it
(usually); otherwise `make avx512-test` runs under Intel SDE. A final job checks that the result hashes
of every function agree across machines and targets.

## License

BSD-3-Clause (see `LICENSE`). The algorithms, constants and scalar routines are derived from Intel's
SVML (Copyright (C) 2023 Intel Corporation, BSD-3-Clause).
