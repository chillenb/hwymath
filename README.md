# hwymath: vectorized double-precision exp

[![CI](https://github.com/chillenb/hwymath/actions/workflows/ci.yml/badge.svg)](https://github.com/chillenb/hwymath/actions/workflows/ci.yml)

Starting from Intel's open-source AVX-512 SVML `exp` (high-accuracy variant), this repository has

1. a bit-identical translation to C with AVX-512 intrinsics, and a bit-identical AVX2 port;
2. faster AVX2 variants with their own polynomials, one of them ("hp") as accurate as Intel's current
   AVX2 HA exp (max 0.511 ulp) at the same speed;
3. a portable [Google Highway](https://github.com/google/highway) version of the hp variant, with run-time
   dispatch, whose results are bit-identical on x86-64 AVX2, AArch64 NEON and SVE/SVE2 (any vector length).

## Getting started

Everything (compilers, Highway, SLEEF, Intel's runtime for comparisons) comes from conda-forge via
[pixi](https://pixi.sh); Intel's assembly is the `third_party/SVML` submodule
([numpy/SVML](https://github.com/numpy/SVML)).

```bash
git clone --recursive https://github.com/chillenb/hwymath && cd hwymath
pixi run make test      # bit-exactness and accuracy checks
pixi run make hwy       # Highway exp on every supported target, with comparisons and benchmarks
pixi run make compare   # x86-64: AVX2 variants vs Intel SVML, glibc, libmvec, SLEEF (accuracy)
pixi run make bench     # x86-64: throughput
pixi run -e cross make xcheck   # AArch64 NEON/SVE/SVE2 under QEMU vs x86-64, bit for bit
```

CI (`.github/workflows/ci.yml`) runs the tests and benchmarks on GitHub's x86-64 and AArch64 runners,
AVX-512 natively when the runner has it (otherwise the AVX-512 correctness tests run under Intel SDE),
and checks that the Highway exp gives bit-identical results on all of them.

Accuracy and speed on an Intel Xeon 6975P-C (AVX-512), a Neoverse-N2 (SVE2) and a Core i5-9500 (AVX2),
compared with Intel SVML, SLEEF, glibc and Highway's own `Exp`: see [Performance.md](Performance.md).

## The SVML translation and AVX2 variants

`svml_exp_d_ha.c` is a translation of `linux/avx512/svml_z0_exp_d_ha.s` from
[numpy/SVML](https://github.com/numpy/SVML) (BSD-3-Clause, Copyright (C) 2023 Intel Corporation) into
C with AVX-512 intrinsics. It keeps the algorithm, the constants (bit-exact, generated from the
assembly into `constants.inc` as hex-float literals), the operation order and every instruction's
rounding mode, so its results should be bit-identical to Intel's.

- `svml_exp8_ha(__m512d)`: the vector routine (needs AVX-512F at run time; only this function is
  compiled for AVX-512 via `__attribute__((target("avx512f")))`).
- `svml_exp_ha_array(x, y, n)`: convenience loop with a masked tail.
- `svml_exp4_ha_avx2(__m256d)`, `svml_exp_ha_array_avx2()` (`svml_exp_d_ha_avx2.c`): AVX2/FMA3 port of
  the same algorithm, designed to be bit-identical to the AVX-512 version. AVX-512-only features are
  replaced by exact equivalents (see the header comment): the round-toward-zero FMA by a
  round-to-nearest FMA plus an exact error-sign correction, `vpermt2pd` by `vpermps` + `blendv` (no
  gathers: they are very slow with the Gather Data Sampling microcode mitigation), `vscalefpd` by
  building 2^M from the bits of t, NaN lanes via the scalar routine.
- `svml_exp4_fast_avx2(__m256d)`, `svml_exp_fast_array_avx2()` (`svml_exp_d_avx2_fast.c`): faster AVX2
  variant that is *not* bit-identical: N rounded to nearest (one FMA, no round-toward-zero emulation),
  2^K-entry tables (default K=3: 8 entries, two `vpermps` + one blend per table), and an own minimax
  polynomial (degree 7 for K=3) fitted by `tools/fit_exp_k3.py` with 60-digit `decimal` arithmetic.
  Special lanes use the same scalar routine. `make FAST_K=2 FAST_DEG=8` builds the 4-entry version.
- `svml_exp4_hp_avx2(__m256d)`, `svml_exp_hp_array_avx2()` (`svml_exp_d_avx2_hp.c`): extra-precision
  variant aiming at ~0.5 ulp. On top of the fast variant: an exact two-part reduction (Ln2Hi with 38
  significant bits, r + r_c carried as a double-double), exp(r) - 1 = r + r^2 Q(r) with Q fitted by
  `tools/fit_exp_k3.py K DEG hp`, and an error-free reconstruction of Th + Th*r (s = fma(Th, r, Th),
  e = fma(Th, r, -(s - Th))), so only the final addition rounds significantly. Default 4-entry
  table (one `vpermps` per lookup, no blends) with Q up to r^9; `make HP_K=3 HP_DEG=8` builds the
  8-entry version (more accurate, somewhat slower; see Performance.md).
- In both non-identical variants the rarely used scalar callout is an out-of-line `cold` function, so
  the hot path needs no stack frame (worth ~10 % per call).
- `svml_dexp_ha_cout_rare()`: the scalar routine for special lanes (overflow, underflow/denormal
  results, ±inf), plain C. Must be compiled without FMA contraction (handled by a pragma).

The header comment of `svml_exp_d_ha.c` explains the algorithm.

## Tests (`make test`)

| Test | Compares with | Status |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dexp_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,002,324 inputs, 0 mismatches |
| `test_exp8` | Intel's `__svml_exp8_ha`: all 8 lanes bit-exact, for both the AVX-512 translation and the AVX2 port | needs AVX-512F (CI: natively or under Intel SDE) |
| `test_model` | scalar model of the AVX-512 instruction semantics (real `fesetround` RZ, `ldexp`, array lookups) vs the AVX2 port | 20,002,324 inputs, 0 mismatches |
| `make insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector path | identical |

Inputs: special values, ±50 ulps around every threshold of both algorithms, table boundaries, and
random values (uniform in [-750, 750] and random bit patterns). The assembly is taken read-only from
`SVML_DIR` (default: the `third_party/SVML` submodule); a copy with the scalar routine made global is assembled into `build/`.

Differences from the assembly that do not affect results: the bitwise ANDs use `vpandq` (AVX-512F)
instead of `vandpd` (AVX-512DQ), and the compiler chooses registers and FMA forms (132/213/231).

## Comparisons

`make compare` and `make bench` (x86-64) and `make hwy` measure these variants against Intel's libsvml,
glibc (scalar and libmvec) and SLEEF; results are in [Performance.md](Performance.md).

`libsvml.so` comes from conda-forge's `intel-cmplr-lib-rt` (closed source, redistributable; used for
comparison only). SLEEF and the Intel runtime are optional for the Makefile: without them they are left
out of the comparisons.

## Highway port (`hwyexp/`)

`hwyexp/hwy_exp.cc` is the hp variant written with [Google Highway](https://github.com/google/highway)
(1.4.0 from conda-forge), with run-time dispatch to the best target:

- targets with true FMA (`HWY_NATIVE_FMA`): the hp algorithm. 4-entry table (one vector on AVX2, two on
  2-lane NEON/SVE-128 via `TwoTablesLookupLanes`, a repeated table for longer vectors; the choice is made
  at run time, so it works for any SVE/RVV vector length), 8-entry table on AVX-512 (one `vpermpd`);
- targets without FMA (SSE2/SSSE3/SSE4): `std::exp`, because the error-free steps would not be exact.
  (Highway's own `Exp` would be ~2x faster there, at ~1 ulp.)

Its results are bit-identical on x86-64 AVX2 and AArch64 NEON/SVE/SVE2 (checked on real hardware in CI,
and locally under QEMU at SVE vector lengths of 128 to 2048 bits with `make xcheck`); the AVX-512
targets, which use the 8-entry table, agree with each other. Accuracy and speed per target:
[Performance.md](Performance.md).

Notes:
- Highway's static AVX2 target also requires AES/PCLMUL, so `-march=x86-64-v3` selects SSSE3 for static
  builds; the cross-check uses `-march=native` for x86-64.
- With a single (static) target, `HWY_DYNAMIC_DISPATCH(F)` expands to `N_<TARGET>::F`: call it inside the
  namespace with an unqualified name.
- Padded lookup tables must be compile-time constants: function-local statics with run-time
  initialization put thread-safe init guards (and calls) into the hot loop, which cost ~25 %.
- The SVML constant macros (`vC1`, ...) collide with identifiers inside Highway's headers, so the shared
  header turns the needed ones into `constexpr` and `#undef`s all of them.

The `cross` pixi environment adds the AArch64 cross compiler and QEMU (`gxx_linux-aarch64`,
`qemu-execve-aarch64`) for `make xcheck`.

## License

BSD-3-Clause (see `LICENSE`). The algorithms, constants and the scalar routine are derived from Intel's
SVML (Copyright (C) 2023 Intel Corporation, BSD-3-Clause).
