# cbrt (double precision)

Source studied: [`svml_z0_cbrt_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_cbrt_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md).

| path | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_cbrt_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py`) |
| `hwy/` | portable Google Highway version of the same algorithm, bit-identical to Intel's on every target with FMA |
| `rcp14_model.h` | exact model of `vrndscalepd(vrcp14pd(m), 0x48)` (eight thresholds), used by the Highway version |
| `tests/` | tests, comparison implementations (Intel SVML, SLEEF), test inputs |
| `tools/` | `gen_svml_constants.py` |
| `cbrt_d.h` | declarations of the C translation |
| `function.mk` | build rules and the `cbrt-*` make targets |

Why it matters: glibc's scalar `cbrt` is up to 3.4 ulp off, and more than half of its results are not
correctly rounded (Performance.md). NumPy's float64 `cbrt` uses it on CPUs without AVX-512.

## Algorithm

With x = ±2^E · m, m ∈ [1, 2) (`vgetexppd`, `vgetmantpd`; both handle denormals):

- **Reciprocal**: rcp ≈ 1/m from `vrcp14pd`, rounded to 4 fractional bits by `vrndscalepd` (imm 0x48), so
  rcp = k/16 with k ∈ {8, …, 16}. R = m·rcp − 1 is computed by one FMA (|R| ≲ 2^-5), and cbrt(m) =
  cbrt(1/rcp) · cbrt(1 + R). Bits 49–52 of rcp (`vpsrlq`) are k − 8, the index into 9-entry tables
  Th + Tl = cbrt(16/k) (16-entry `vpermt2pd` tables, double-double).
- **Exponent**: K = ⌊E/3⌋ and j = E − 3K. ES = E + 1.5·2^52 (an integer in the low mantissa bits),
  t = ES·(1/3)⁺ − 2^51 = E/3 + 1/4 + tiny (one FMA; 1/3 is rounded up), K = ⌊t⌋ (`vrndscalepd` imm 9),
  and ES − 3K = 1.5·2^52 + j (FMA). The low bits of that are the index into Sh + Sl = 2^(j/3)
  (`vpermpd`). The 1/4 offset keeps E/3's fractional parts {0, 1/3, 2/3} away from integers, so ⌊t⌋ is
  exact for every E.
- **Result**: ShTh = RN(Sh·Th). Low = (Sh·Th − ShTh) + Sh·Tl + Sl·Th (three FMAs). P(R) has degree 9
  with cbrt(1 + R) ≈ 1 + R·P(R). Then cbrt(|x|) = 2^K · (ShTh + (Low + (ShTh·R)·P(R))), with the 2^K
  from `vscalefpd`, and the sign of x is ORed in (`vorpd`).
- **No callout.** 0 (E = −∞ → 2^−∞ scaling), ±inf and NaN come out right by themselves, so the vector
  path has no branch. The scalar routine `__svml_dcbrt_ha_cout_rare_internal` is VML's scalar fallback.
  It uses a 32-entry table rc ≈ m^(−1/3) (8 bits) and cbrt(m) = m·rc²·(1 − d)^(−2/3), d = 1 − m·rc³, in
  double-double arithmetic with Dekker splits. It handles denormals by scaling with 2^300 and integer
  division of the exponent by 3.

Details: the header comment of `svml/svml_cbrt_d_ha.c`.

## The C translation

`svml/svml_cbrt_d_ha.c` translates the assembly into C with AVX-512 intrinsics. It keeps the algorithm, the
constants (bit-exact hex-float literals), the operation order and every instruction's rounding mode. Its
results are bit-identical to Intel's (`test_cbrt8`).

- `svml_cbrt8_ha(__m512d)`: the vector routine (needs AVX-512F at run time; only it and the array loop
  are compiled with `__attribute__((target("avx512f")))`). `svml_cbrt_ha_array(x, y, n)`: loop with a
  masked tail.
- `svml_dcbrt_ha_cout_rare()`: the scalar routine, plain C, compiled without FMA contraction. It is about
  0.50 ulp but slow (36 ns/value); the vector routine never calls it.
- Differences from the assembly that do not affect results: the bitwise operations are `vpandq`/`vporq`
  (AVX-512F) instead of `vandpd`/`vorpd` (AVX-512DQ). The tables are read through an empty `asm`
  barrier, because GCC otherwise replaces the load of `vTl`'s all-zero upper half with a zeroing
  instruction. The compiler also chooses registers and FMA forms (132/213/231).

## Highway version (`hwy/`)

`hwy/hwy_cbrt.cc` ports the same algorithm with run-time dispatch. On targets with true FMA
(`HWY_NATIVE_FMA`), each AVX-512-only step is replaced by an exact equivalent, so the results are
**bit-identical to Intel's `__svml_cbrt8_ha`** (as executed by Intel CPUs / Intel SDE). Targets without
FMA (SSE2/SSSE3/SSE4, EMU128, SCALAR) fall back to `std::cbrt`.

| AVX-512 feature | emulation |
|---|---|
| `vgetexppd`, `vgetmantpd` | exponent and mantissa bits. Denormals are first multiplied by 2^64 (exact) and the exponent corrected. ES = 1.5·2^52 + E is built directly as bits (`vaddpd` is exact there) |
| `vrcp14pd` + `vrndscalepd` imm 0x48 | k = 16 − #{i : m ≥ T_i}, eight compares against the thresholds of `rcp14_model.h`; rcp = k/16 exactly (see below) |
| `vrndscalepd` imm 9 (floor) | `Floor` |
| `vpermpd` (2^(j/3), 3 of 8 entries used) | `TableLookupLanes` on a 4-entry table (two-table lookup on 2-lane vectors) |
| `vpermt2pd` (cbrt(16/k), 9 of 16 entries used) | `TableLookupLanes` (≥ 8 lanes) / `TwoTablesLookupLanes` (4 lanes) / two 4-entry lookups + select (2 or 6 lanes) on entries k = 8…15, and a select for k = 16 (Th = 1, Tl = 0) |
| `vpsrlq` index from rcp's bits | k − 8 = low 3 bits of the counter (k = 16 is the select) |
| `vscalefpd` | multiplication by 2^K built from exponent bits: K ∈ [−358, 341] and the scaled value is in [1, 2), so the result is always normal and the product is exact |
| FMAs/multiplies with `{rn-sae}` | `MulAdd`/`MulSub`/`NegMulAdd`/`Mul` (round to nearest by default) |
| 0, ±inf, NaN (handled implicitly by the instructions) | `x + x`, which is exactly what the AVX-512 code returns, NaN payloads included (checked under SDE) |

**`vrcp14pd`.** Its output is implementation-defined (documented relative error < 2^-14), so the question
was whether rounding to 4 fractional bits makes the result independent of the approximation. **It does
not**: rcp14 ≈ 1/m is only accurate to about 2^-14, and the rounding of 16·rcp14 to an integer is
ambiguous wherever 16/m is within ~2^-10 of a half-integer, i.e. near m = 32/(2k+1). Probing Intel's
`vrcp14pd` under SDE showed three things:

- For m ∈ [1, 2) the output depends only on the **top 16 mantissa bits** of m (except rcp14(1) = 1
  exactly) and has 16 fractional bits.
- It is monotone. So the rounded value is a step function of the 16-bit prefix, switching at eight
  thresholds T_i, all multiples of 2^-16: 0x1.0842p0, 0x1.1a7dp0, 0x1.2f69p0, 0x1.47adp0, 0x1.642cp0,
  0x1.8619p0, 0x1.af29p0, 0x1.e1e2p0.
- Exact rounding of 16/m would switch at 32/(2k+1) instead, which differs from these by −1.1 to +1.4
  prefixes. That is about 1 mantissa in 12,850 with a different table index. About 1.1 % of those, i.e.
  ~1 in 10^6 random inputs, would give a different final result (either closer to or farther from the
  true root).

The model `rcp14_model.h` (k = 16 − #{i : |m| ≥ T_i}) is therefore exact, and `test_rcp14` checks it
**exhaustively**: all 2^16 prefixes × both signs × low bits zero / all ones / one / random, plus 4M random
mantissas, against `vrndscalepd(vrcp14pd(m), 0x48)` under SDE. It also checks that rcp14 depends on no
more than the 16 prefix bits.

The thresholds were measured on Intel SDE's `vrcp14pd`. `test_rcp14` has since passed natively in CI on an
Intel Xeon 6975P-C, and the same instruction follows Intel's reference algorithm exactly on AMD Zen 4 and Zen 5
(EPYC 9V74, 9V45: `test_approx14`, `verify_rcp14`), so the thresholds hold there too. On AVX-512 targets the
Highway version now uses the instruction itself (`Round(16·vrcp14pd(m))/16`); the thresholds serve the
other targets.

Performance notes:
- A fast path for vectors whose lanes are all normal skips the denormal scaling and the final special-value
  select (2.75 → 2.28 ns/value on AVX2). Other vectors take the general path. Neither path needs the
  scalar routine.
- Padded lookup tables are `constexpr` (the generated `.inc` tables become `constexpr` in C++ via
  `CBRT_TABLE_CONST`), so the hot loop has no run-time-initialized statics.
- No gathers (`GatherIndex` is only reachable with 1-lane vectors, which no x86 FMA target has).

## Tests (`make cbrt-test`, `make cbrt-avx512-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dcbrt_ha_cout_rare_internal` (made global with `sed`): results, return codes, FP exception flags | 20,019,282 inputs, 0 mismatches |
| `test_cbrt8` | Intel's `__svml_cbrt8_ha`, all 8 lanes bit-exact | 2,019,288 inputs under Intel SDE (`-spr`), 0 mismatches; needs AVX-512F |
| `test_rcp14` | `vrndscalepd(vrcp14pd(m), 0x48)` vs `rcp14_model.h`: all 2^16 mantissa prefixes and random mantissas | 4,524,288 inputs under SDE, 0 mismatches; rcp14 depends on the 16 prefix bits only |
| `hwy_test` | every Highway target: max error < 1 ulp, special values, odd symmetry bit for bit, agreement between targets, and with AVX-512F (SDE) bit-identity with the C translation | AVX3_SPR, AVX3_DL, AVX3, AVX2: 0 of 2,019,282 results differ from the translation (SDE) |
| `make cbrt-insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector routine vs the assembly | identical |

Inputs (`tests/inputs.h`):
- special values (±0, ±inf, quiet and signalling NaNs with various payloads);
- every power of two from 2^-1074 to 2^1023 with its neighbours, both signs;
- perfect cubes n³ (n ≤ 1536) scaled by 2^(3s) into the normal, tiny, denormal and huge ranges;
- ±20 ulps around each rcp14 threshold and each exact boundary 32/(2k+1), with exponents covering
  E mod 3 = 0, 1, 2 and denormals;
- ±5 ulps around the scalar routine's table boundaries 1 + i/32;
- random values: arbitrary bit patterns, uniform mantissas with small exponents, uniform in
  [−1000, 1000].

`hwy_test` also reports that all perfect cubes and powers of 8 among the inputs get exact cube roots.
`make cbrt-hashes` / `cbrt-avx512-test` write `arch cbrt target svml n hash` lines, plus the C
translation's own hash (`svml-c-avx512`) when AVX-512F is available. Every `svml` line must agree across
targets and machines (`common/tools/check_hashes.py`). That extends the bit-identity check to
AArch64/SVE via CI. A one-off QEMU cross-check gave the translation's hash for the 2,019,282 test
inputs on AArch64 NEON and on SVE and SVE2 at 128-, 256-, 384- (6 lanes: split 8-entry lookup), 512- and
2048-bit vectors (static builds with `-DHWY_COMPILE_ONLY_STATIC`, like `exp`'s `tools/xcheck.sh`; not a
make target here).

## Comparisons

`make cbrt-hwy` measures the Highway version against Intel's libsvml (conda-forge `intel-cmplr-lib-rt`,
closed source, used for comparison only: `__svml_cbrt4_ha_l9`, `__svml_cbrt8_ha_z0`), SLEEF u10, the C
scalar routine and glibc's scalar `cbrt` (glibc 2.28's libmvec has no cbrt). Results are in
[Performance.md](Performance.md).
