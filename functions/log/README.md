# log (double precision)

Source studied: [`svml_z0_log_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_log_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md).

This directory also holds the parts shared by the whole family (log, [log2](../log2/README.md),
[log10](../log10/README.md), [log1p](../log1p/README.md)): the four SVML routines use the same reduction
(rounded `vrcp14pd`, 16-entry tables), so the emulation of the AVX-512 operations, the analysis of
`vrcp14pd`, the constant generator and the test driver exist once, here. Each function still has its own
directory, because the C translations, callouts, polynomials and inputs differ (log2's callout is an
entirely different algorithm), and the per-function make targets stay independent.

| directory | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_log_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py log`) |
| `hwy/` | Google Highway version (`hwy_log.cc`); **shared**: `log_hwy-inl.h` (emulation of the AVX-512 operations), `rcp14_thresholds.h` (the rounded-`vrcp14pd` model), `log_tables.h` |
| `tests/` | tests and test inputs; **shared**: `log_family_test.h` (Highway test/benchmark driver), `log_family_xcheck.h`, `extras*` (Intel SVML, SLEEF, libmvec) |
| `tools/` | **shared**: `gen_svml_constants.py FUNCTION [--cxx]`, `verify_rcp14.c` (the `vrcp14pd` model against the reference and the CPU or Intel SDE), `xcheck.sh FUNCTION` (QEMU) |
| `log_d.h` | declarations of the C translation |
| `function.mk` | build rules and the `log-*` make targets |

## Algorithm

Vector path (8 lanes): x = 2^E · m with m ∈ [1, 2) by `vgetexppd`/`vgetmantpd` (exact, also for
denormals). `rcp = vrndscalepd(vrcp14pd(m))` rounds the approximate reciprocal to a multiple of 1/32:
rcp = k/32, k ∈ 16..32. Then R = rcp·m − 1 is exact (one FMA, |R| ≲ 2^-5), and
log(x) = E·ln2 − log(rcp) + log(1 + R). For rcp < 0.75, E is incremented and the table holds −log(2·rcp),
so the table index is just the top 4 mantissa bits of rcp (`vpsrlq` 48, low 4 bits used by `vpermt2pd`;
k = 32 and k = 16 share entry 0 = 0). Th + Tl = −log(rcp) as a double-double (16-entry tables). The sum is
evaluated as High = E·L2H + Th (exact, L2H has 37 bits), HR = High + R with its rounding error Rl
(Fast2Sum), and log(1+R) − R = R²·(c2 + c3 R + … + c9 R⁷) in an Estrin-like scheme:
result = HR + ((R²·p + Rl) + (E·L2L + Tl)).

Lanes with x ≤ 0, ±inf or NaN (`vfpclasspd` 0xdf; denormals stay in the vector path) are recomputed by the
scalar routine `__svml_dlog_ha_cout_rare_internal`, which also raises the FP flags and returns an error code
(1: x < 0 or −inf → NaN, invalid; 2: ±0 → −inf, divide-by-zero). Its own algorithm: denormals scaled by
2^60, |x − 1| ≤ 2^-6.68 by a degree-8 polynomial, otherwise a 65-entry table {rcp, −log(rcp) hi, lo} indexed
by the mantissa rounded to 6 bits, with the reduced argument split (y rounded to 31 bits) so that
rcp·yh − 1 is exact. Details: the header comment of `svml/svml_log_d_ha.c`.

## The C translation

`svml/svml_log_d_ha.c` translates the assembly into C with AVX-512 intrinsics, keeping the algorithm, the
constants (bit-exact hex-float literals), the operation order and every instruction's rounding mode, so its
results are bit-identical to Intel's (`test_log8`), and the instruction mix of the vector path is identical
(`make log-insns`).

- `svml_log8_ha(__m512d)`: the vector routine (AVX-512F + DQ for `vfpclasspd`/`vxorpd`; only this function
  is compiled for AVX-512 via `__attribute__((target(...)))`). `svml_log_ha_array(x, y, n)`: loop with a
  masked tail (padding lanes are 1.0, not 0, so that they never take the callout).
- `svml_dlog_ha_cout_rare()`: the scalar routine, plain C (no FMA contraction). The NaNs it creates
  (0/0) are returned with the sign of x86's default NaN on every architecture (`copysign(z/z, -1.0)`; on
  x86 a no-op), so that results are identical on AArch64, whose default NaN is positive.

## Highway version (`hwy/`)

`hwy/hwy_log.cc` is the same algorithm written with [Google Highway](https://github.com/google/highway)
(1.4.0 from conda-forge), with run-time dispatch. On every target with true FMA (`HWY_NATIVE_FMA`) it
performs the same floating-point operations in the same order with the same (round-to-nearest) rounding
as `__svml_log8_ha`, and its results are **bit-identical** to Intel's AVX-512 code: verified under Intel
SDE for AVX2, AVX3, AVX3_DL, AVX3_SPR (`-spr`) and AVX10_2 (`-dmr`) against the C translation, and under
QEMU for AArch64 NEON, SVE and SVE2 at 128/256/512/2048-bit vectors against x86-64 (`make log-xcheck`).
Targets without FMA (SSE2/SSSE3/SSE4) use `std::log`.

The AVX-512-only operations, emulated exactly on every target except the AVX-512 ones (`hwy/log_hwy-inl.h`):

| AVX-512 | emulation |
|---|---|
| `vgetmantpd` (interval [1,2)), `vgetexppd` | mantissa bits OR 1.0; exponent field OR'ed into the bits of 2^52 and 2^52 + 1023 subtracted (exact, no int→double conversion). Denormals: pre-scaled by 2^52, only in a slower variant of the vector code used for vectors that contain a denormal or special lane |
| `vrcp14pd` + `vrndscalepd` (to 1/32) | exact model of the rounded reciprocal: 16 thresholds, a 32-entry `uint32` table and integer operations on the high word of m (below) |
| `vpermt2pd` (16-entry tables) | `TwoTablesLookupLanes` (one `vpermt2pd` on AVX-512); 4 lanes: two of them plus a select; 2 lanes: four plus three selects; ≥ 16 lanes: `TableLookupLanes` on a table repeated at compile time. No gathers |
| `vfpclasspd` 0xdf + callout | `!(x >= DBL_MIN && x < inf)` decides between the fast vector code and the slow one followed by the scalar callout for the lanes with `!(x > 0 && x < inf)` |
| masked `vaddpd` (E + 1 if rcp < 0.75) | `MaskedAddOr` |
| masked `vxorpd` (Rh = 0 in special lanes) | not needed: those lanes are overwritten by the callout |
| `{rn-sae}` / `{sae}` | all rounding is round-to-nearest, so plain `MulAdd`/`MulSub`/`Add`/`Mul` (fused on FMA targets) |

On the AVX-512 targets (AVX3, AVX3_DL, AVX3_ZEN4, AVX3_SPR, AVX10_2; `HWYMATH_NATIVE_AVX512` from
`common/include/hwymath/hwy/avx512_native-inl.h`) the helpers use Intel's instructions instead: `vgetexppd` and
`vgetmantpd $8` {sae} (exact for denormals too, so there is only one vector path), `vrcp14pd` +
`vrndscalepd $0x58` {sae}, the table index `vpsrlq $48` of rcp without masking (as in Intel's code: `vpermt2pd`
only reads the low 4 bits; `Lookup16` calls `_mm512_permutex2var_pd` directly so that the unmasked index is
legitimate) and `vfpclasspd $0xdf` for the special-lane test. The AVX3 loop then has 29 vector/mask instructions
besides loads, stores and register moves (37 with the emulations), as many as Intel's main path less its
`vxorpd` (30, counting `kmovw` like our `kortestb`).

### `vrcp14pd`: why exact emulation is possible

`vrcp14pd` returns an approximation of 1/m with relative error < 2^-14; the SDM leaves its bits to the
implementation, but Intel publishes reference C code (`recip14.c`, RCP14D: a 64-segment piecewise-linear
formula on the top 16 fraction bits, truncated to 17 bits; not the correctly rounded value — about 25 % of
results equal 1/m rounded to 16 bits, errors up to 5.6 units of 2^-17). `common/include/hwymath/rcp14_ref.h` is an integer
rewrite of it (copied from the tan port, branch `port/trig`). SVML rounds the result at once to 5
fractional bits (`vrndscalepd` imm 0x58: nearest, ties to even, scale 2^-5), so only the 16 places where the
rounded value changes matter:

- k(m) = 32·rcp is a non-increasing step function of m with exactly 16 steps, at the mantissas
  `kRcpThreshold[]` of `hwy/rcp14_thresholds.h`. They are multiples of 2^-16 (the reference depends only on
  the top 16 fraction bits, and on m = 1 exactly) and lie up to ~2 units of 2^-16 away from the exact
  boundaries 64/(2k+1), on either side — so RN(32/m) with a division would **not** be exact.
- Proof, given the reference: `tools/verify_rcp14.c` rounds the reference for all 2^16 fraction prefixes
  (which covers every m) and compares with the step function; it also checks the emulation (below) for all
  2^20 high words of m, i.e. for every m. Both run on any machine (`make log-test`).
- The reference against the CPU (`make log-avx512-test` / `log-rcp14-check`, natively on AVX-512 hardware,
  e.g. in CI, or under Intel SDE): `vrcp14pd` equals the reference bit for bit for all 2^16 prefixes with 8
  low-bit patterns each, and k(m) from the CPU matches the step function at the thresholds found again by
  bisection, at every double near them, at both ends of all 2^24 intervals covering [1, 2) and at 16 M
  random mantissas. Under SDE 10.13 `-spr`: 0 mismatches, also in a longer run (`verify_rcp14 22 36 13`:
  134 M doubles around the thresholds and 268 M interval ends within 2^-13 of them, 2.5 min); 16 M random
  full mantissas also matched the reference bit for bit.
- The emulation: each bin of m with the same top 5 mantissa bits (width 1/32) contains at most one threshold
  (they are ≥ 0.0333 apart), so `kRcpBins[b]` = high word of the bin's threshold | (k above it − 16), and
  k − 16 = (entry & 15) + ((hi32(m) | 15) < entry) — a 32-entry `uint32` lookup (one `vpermt2d` on
  AVX-512) and a few integer operations; rcp's bits are 0x3fe0000000000000 + ((k − 16) << 48) (k = 32
  gives 1.0) and the table index is (k − 16) & 15. This is cheaper than evaluating the full reference
  formula and rounding it.

So the Highway version reproduces Intel's reference (and SDE's) `vrcp14pd`. If another AVX-512 implementation
(e.g. AMD Zen 4/5) rounds differently near these 16 thresholds, Intel's own SVML binary would give different
results there on that CPU, while this version stays identical on every machine (the CI hash cross-check).
An earlier version computed k from a quadratic estimate of 32/m and a 16-entry threshold lookup (also exact,
but 4.7 ns/value instead of 2.5 on AVX2).

Other notes:
- The vector code does not raise FP exceptions like the `{sae}` assembly (Highway's operations may set
  inexact etc.); special lanes go through the same C callout, which raises the documented flags.
- The constants (`hwy/hwy_log_constants.inc`) are generated from the assembly as C++ `constexpr` (`--cxx`),
  so that the padded lookup tables are built at compile time.
- Padding lanes of the tail are zeros: they take the slow path, whose results are discarded.

## Tests (`make log-test`, `make log-avx512-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dlog_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,029,320 inputs, 0 mismatches |
| `test_log8` | Intel's `__svml_log8_ha`: all 8 lanes bit-exact | 2,029,320 inputs, 0 mismatches (Intel SDE `-spr`); needs AVX-512F/DQ |
| `hwy_test` | every Highway target: error < 1 ulp, special values vs `logl`, agreement between targets, special lanes bit-identical to the C callout; with AVX-512: every target bit-identical to the C translation | all checks pass (AVX2 natively; AVX2, AVX3, AVX3_DL, AVX3_SPR under SDE `-spr`, also AVX10_2 under `-dmr`) |
| `make log-insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector path vs the assembly | identical |
| `verify_rcp14` (in `log-test`, `log-avx512-test`, `log-rcp14-check`) | the `vrcp14pd` model and its emulation vs Intel's reference algorithm (any CPU) and vs `vrcp14pd` (AVX-512) | 0 mismatches |
| `make log-xcheck` | Highway on AArch64 NEON/SVE/SVE2 (QEMU, 128–2048-bit SVE) vs x86-64 | bit-identical |

Inputs (`tests/inputs.h`): special values (signed zeros, infinities, NaNs with payloads, denormals, DBL_MAX),
±50 ulps around every threshold of the vector path and the callout (including the log2 callout's windows
around 1, so that the file serves log2 and log10 unchanged), ±8 ulps around the 16 `vrcp14pd` thresholds and
the 129 boundaries 1 + i/128 at 11 exponents (denormal to 2^1023), random values in [0, 4), near 1
(uniform in [0.98, 1.02) and 1 ± 2^-[1,60]), random positive bit patterns and random bit patterns. The
assembly is read from the submodule; a copy with the scalar routine made global is assembled into
`build/log/`. `make log-hashes` writes `arch log <target> svml n hash` lines (and `C-AVX512` for the
translation when the CPU has AVX-512): all of them must be equal across targets and machines.

Differences from the assembly that do not affect results: the compiler chooses registers and FMA forms
(132/213/231), and the order of independent instructions.
