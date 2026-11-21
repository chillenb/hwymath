# exp (double precision)

Source studied: [`svml_z0_exp_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_exp_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md).

| directory | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_exp_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py`) |
| `avx2/` | bit-identical AVX2 port; faster, non-identical "fast" and "hp" variants with their own polynomials |
| `hwy/` | portable Google Highway versions: of the hp variant (`hwy_exp.cc`), and of Intel's algorithm, bit-identical to the assembly (`hwy_exp_svml.cc`) |
| `tests/` | tests, benchmarks, comparison implementations (Intel SVML, SLEEF, libmvec), test inputs |
| `tools/` | `gen_svml_constants.py`, `fit_poly.py` (minimax fits, 60-digit `decimal`), `xcheck.sh` (QEMU) |
| `exp_d.h` | declarations of all the C kernels |
| `function.mk` | build rules and the `exp-*` make targets |

## Algorithm

exp(x) = 2^M · 2^(j/16) · exp(r): N = M + j/16 is x·log2(e) truncated to a multiple of 1/16 by one
round-toward-zero FMA with a shifter (1.5·2^48, whose low bits hold the exponent bias), r = x − N·ln2
by Cody–Waite reduction, 2^(j/16) = Th·(1 + Tlr) from 16-entry tables (`vpermt2pd`), exp(r) − 1 by a
degree-6 polynomial, and the scaling by 2^M with `vscalefpd`. Lanes with |x| ≥ 707.703 (and ±inf) are
recomputed by a scalar routine with its own 64-entry table, which raises the floating-point flags and
handles results in the denormal range with a single rounding. Details: the header comment of
`svml/svml_exp_d_ha.c`.

## The C translation and the AVX2 versions

`svml/svml_exp_d_ha.c` translates the assembly into C with AVX-512 intrinsics. It keeps the algorithm,
the constants (bit-exact hex-float literals), the operation order and every instruction's rounding mode,
so its results are bit-identical to Intel's (`test_exp8`).

- `svml_exp8_ha(__m512d)`: the vector routine (needs AVX-512F at run time; only this function is
  compiled for AVX-512 via `__attribute__((target("avx512f")))`). `svml_exp_ha_array(x, y, n)`: loop
  with a masked tail.
- `svml_dexp_ha_cout_rare()`: the scalar routine for special lanes (overflow, underflow/denormal
  results, ±inf), plain C. Must be compiled without FMA contraction (handled by a pragma).
- `avx2/svml_exp_d_ha_avx2.c`, `svml_exp4_ha_avx2(__m256d)`: AVX2/FMA3 port, bit-identical to the
  AVX-512 version. AVX-512-only features are replaced by exact equivalents (see its header comment): the
  round-toward-zero FMA by a round-to-nearest FMA plus an exact error-sign correction, `vpermt2pd` by
  `vpermps` + `blendv` (no gathers: they are very slow with the Gather Data Sampling microcode
  mitigation), `vscalefpd` by building 2^M from the bits of t, NaN lanes via the scalar routine.
- `avx2/exp_d_avx2_fast.c`, `svml_exp4_fast_avx2(__m256d)`: faster, *not* bit-identical: N rounded to
  nearest (one FMA), 2^K-entry tables (default K=3: two `vpermps` + one blend per table) and an own
  minimax polynomial (degree 7 for K=3). `make EXP_FAST_K=2 EXP_FAST_DEG=8` builds the 4-entry version.
- `avx2/exp_d_avx2_hp.c`, `svml_exp4_hp_avx2(__m256d)`: extra precision, ~0.5 ulp. On top of the fast
  variant: an exact two-part reduction (Ln2Hi with 38 significant bits, r + r_c as a double-double),
  exp(r) − 1 = r + r²·Q(r) with Q fitted by `tools/fit_poly.py K DEG hp`, and an error-free reconstruction
  of Th + Th·r (s = fma(Th, r, Th), e = fma(Th, r, −(s − Th))), so only the final addition rounds
  significantly. Default 4-entry table (one `vpermps` per lookup) with Q up to r⁹;
  `make EXP_HP_K=3 EXP_HP_DEG=8` builds the 8-entry version (more accurate, somewhat slower).
- In both non-identical variants the rarely used scalar callout is an out-of-line `cold` function, so
  the hot path needs no stack frame (worth ~10 % per call).

## Highway versions (`hwy/`)

### hp variant

`hwy/hwy_exp.cc` is the hp variant written with [Google Highway](https://github.com/google/highway)
(1.4.0 from conda-forge), with run-time dispatch to the best target:

- targets with true FMA (`HWY_NATIVE_FMA`): the hp algorithm. 4-entry table (one vector on AVX2, two on
  2-lane NEON/SVE-128 via `TwoTablesLookupLanes`, a repeated table for longer vectors; chosen at run
  time, so it works for any SVE/RVV vector length), 8-entry table on AVX-512 (one `vpermpd`);
- targets without FMA (SSE2/SSSE3/SSE4): `std::exp`, because the error-free steps would not be exact
  (Highway's own `Exp` would be ~2x faster there, at ~1 ulp).

Its results are bit-identical on x86-64 AVX2 and AArch64 NEON/SVE/SVE2 (on real hardware in CI, and
under QEMU at SVE vector lengths of 128 to 2048 bits with `make exp-xcheck`); the AVX-512 targets, which
use the 8-entry table, agree with each other.

Notes:
- Highway's static AVX2 target also requires AES/PCLMUL, so `-march=x86-64-v3` selects SSSE3 for static
  builds; the cross-check uses `-march=native` for x86-64.
- With a single (static) target, `HWY_DYNAMIC_DISPATCH(F)` expands to `N_<TARGET>::F`: call it inside the
  namespace with an unqualified name.
- Padded lookup tables must be compile-time constants: function-local statics with run-time
  initialization put thread-safe init guards (and calls) into the hot loop, which cost ~25 %.
- The SVML constant macros (`vC1`, ...) collide with identifiers inside Highway's headers, so
  `hwy/hwy_exp_constants.h` turns the needed ones into `constexpr` and `#undef`s all of them.

### Intel's algorithm (`hwy_exp_svml.cc`)

`hwy/hwy_exp_svml.cc` (`hwy_exp_svml_array`) ports the algorithm of the assembly itself, like
`functions/exp2/hwy/`. On every target with true FMA (`HWY_NATIVE_FMA`: AVX3*, AVX2, NEON, SVE, SVE2 at
any vector length) its results are bit-identical to Intel's `__svml_exp8_ha`; targets without FMA
(SSE2/SSSE3/SSE4) use `std::exp`. The AVX-512-only operations are emulated exactly by small helper
functions (to be replaced by the native instructions on AVX-512 targets):

- `{rz-sae}` FMA (`ShifterRZ`): as in the AVX2 port, t = RN(x·log2(e) + S), the sign of its exact
  rounding error from a second FMA, and one ulp down on the bit pattern where RN rounded up (a masked
  `vpsubq` on AVX-512, an add of the all-ones compare mask elsewhere);
- `vpermt2pd` (`Lookup16`): `vpermt2pd` itself on AVX-512 (no index masking: it reads the low 4 bits of
  t); on AVX2 two `TwoTablesLookupLanes` (`vpermd` + `blendv`) for entries 0–7 and 8–15 and a `blendv`
  on bit 3 of the index shifted into the sign bit; selects of two-table lookups on 2-lane vectors; a
  repeated table and `TableLookupLanes` for ≥ 16 lanes. No gathers;
- `vscalefpd` (`Scalef`): bits 4–15 of t are 1023 + floor(N), so 2^floor(N) is two shifts away; it is
  normal and the product exact for every lane that the scalar routine does not redo;
- Intel's `and` of r with `0xbfff…` (it keeps |r| < 2 for special inputs) is a no-op for every lane whose
  vector result is kept, so it is left out;
- lanes with |x| ≥ 707.703 *or NaN* go to the C scalar routine `svml_dexp_ha_cout_rare` (out of line,
  `cold`); for NaN it returns x·x, the quieted NaN that the AVX-512 code produces.

Checked by `hwy_test` (bit-identity of every svml target with the C AVX2 port, and with AVX-512F —
natively or under Intel SDE — with the C/AVX-512 translation and Intel's assembly), by `make exp-xcheck`
(QEMU) and by the CI hash cross-check (variant `svml`; the C AVX2 port, the AVX-512 translation and the
assembly write hash lines of the same variant).

## Tests (`make exp-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dexp_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,002,324 inputs, 0 mismatches |
| `test_exp8` | Intel's `__svml_exp8_ha`: all 8 lanes bit-exact, for the AVX-512 translation and the AVX2 port | 100,002,328 inputs, 0 mismatches (CI, Xeon 6975P-C); needs AVX-512 (else Intel SDE) |
| `test_model` | scalar model of the AVX-512 instruction semantics (real `fesetround` RZ, `ldexp`, array lookups) vs the AVX2 port | 20,002,324 inputs, 0 mismatches |
| `hwy_test` | every Highway target, both versions: error bounds, agreement between targets, special values; hp: bit-identity with the C hp variant; svml: bit-identity with the C AVX2 port, and with AVX-512F with the C/AVX-512 translation and `__svml_exp8_ha` | all checks pass (AVX3_SPR/AVX3_DL/AVX3/AVX2 under SDE `-spr`) |
| `make exp-xcheck` | both Highway versions as static AArch64 NEON/SVE/SVE2 binaries under QEMU (SVE at 128–2048 bits) vs x86-64 AVX2 | bit-identical |
| `make exp-insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector path vs the assembly | identical |

Inputs (`tests/inputs.h`): special values, ±50 ulps around every threshold of both algorithms, table
boundaries, and random values (uniform in [−750, 750] and random bit patterns). The assembly is read
from the submodule; a copy with the scalar routine made global is assembled into `build/exp/`.

Differences from the assembly that do not affect results: the bitwise ANDs use `vpandq` (AVX-512F)
instead of `vandpd` (AVX-512DQ), and the compiler chooses registers and FMA forms (132/213/231).

## Comparisons

`make exp-hwy`, and on x86-64 `make exp-compare` / `make exp-bench`, measure these versions against
Intel's libsvml (conda-forge `intel-cmplr-lib-rt`, closed source, used for comparison only), glibc
(scalar and libmvec), SLEEF and Highway's own `Exp`; results in [Performance.md](Performance.md).
