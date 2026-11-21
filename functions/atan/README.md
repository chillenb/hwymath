# atan (double precision)

Source studied: [`svml_z0_atan_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_atan_d_ha.s)
(AVX-512, high accuracy; entry `__svml_atan8_ha`). Accuracy and speed: [Performance.md](Performance.md).
The scalar routine is translated mechanically like asin's ([../asin/README.md](../asin/README.md)).

| path | contents |
|---|---|
| `svml/svml_atan_d_ha.c` | bit-identical C/AVX-512 translation of the vector routine |
| `svml/svml_atan_d_ha_constants.inc` | vector-path constants and tables, generated from the assembly by `tools/gen_svml_constants.py` |
| `svml/svml_datan_ha_cout_rare.inc` | the scalar routine, translated mechanically by `functions/asin/tools/asm2c.py` |
| `hwy/` | portable Google Highway version |
| `tests/` | tests, benchmark, comparison implementations (Intel SVML, SLEEF), test inputs |
| `atan_d.h` | declarations of the C kernels |
| `function.mk` | build rules and the `atan-*` make targets |

## Algorithm

With a = |x| (the sign of x is restored with a final xor):

- a < 7.875: t = RN(a + 1.5·2⁵⁰) rounds a to a multiple of 1/4 (ties to even); j = 0…31 are the low
  mantissa bits of t and x₀ = j/4 = t − 1.5·2⁵⁰. Diff = a − x₀ (`vreducepd`, exact, |Diff| ≤ 1/8),
  Y = 1 + x₀·a (one FMA), and atan(a) = atan(j/4) + atan(Diff/Y), with atan(j/4) = AtanHi[j] + AtanLo[j]
  from two 32-entry tables, each read with two `vpermt2pd` (entries 0–15 and 16–31) and a blend on j ≥ 16;
- a ≥ 7.875 (and ±inf, NaN): Diff = −1, Y = min(a, 2¹²⁸), atan(a) = π/2 + atan(−1/Y).

The quotient is a double-double R + R_low: Rcp = `vrcp14pd`(Y) (relative error < 2⁻¹⁴),
eps = 1 − Y·Rcp, Rcp₂ = Rcp(1 + eps)(1 + eps²) (two FMAs), R = RN(Diff·Rcp₂),
R_low = (Diff·Rcp₂ − R) + R·(1 − Y·Rcp₂), minus (x₀a − (Y − 1))·Rcp₂·R for a < 7.875 (the rounding error
of Y). Finally atan(a) = (AtanHi + R) + (R³·P(R²) + (R_low + AtanLo + err)), where err is the rounding
error of AtanHi + R and P has degree 5. There is **no scalar callout** in the vector path: ±inf gives
±RN(π/2), NaN the quieted input. The assembly still contains the scalar routine
`__svml_datan_ha_cout_rare_internal` (a complete table-driven scalar atan, return code always 0); it is
translated and tested on its own (`test_rare`) but not used.

## The C translation

`svml/svml_atan_d_ha.c` keeps the algorithm, the constants, the operation order and every instruction's
rounding mode: `svml_atan8_ha(__m512d)` (bit-identical to Intel's, `test_atan8`; needs AVX-512F and DQ for
`vreducepd`), `svml_atan_ha_array` and `svml_datan_ha_cout_rare(&x, &r)` (the generated scalar
routine, portable C). As in asin, empty `asm` barriers keep GCC from rewriting the sign operations as
`and`/`vpternlogq`, so that `make atan-insns` compares equal instruction mixes.

## Highway version (`hwy/`)

| AVX-512 feature | emulation |
|---|---|
| `{rn-sae}` arithmetic and FMAs | the same IEEE operations (`Mul`, `Add`, `Sub`, `MulAdd`, `MulSub`, `NegMulAdd`) |
| `vreducepd $0x28` (a − round(4a)/4, ties to even) | a − (RN(a + 1.5·2⁵⁰) − 1.5·2⁵⁰): the same rounding, and the subtraction is exact (Sterbenz, or x₀ = 0) |
| `vpermt2pd` on the 32-entry tables (4 per call) | `TableLookupLanes` on a repeated table when a vector holds 32+ doubles, else `TwoTablesLookupLanes` on chunks of 2N entries merged with `IfThenElse` on the high index bits (AVX-512: 2 per table, as Intel; AVX2: 4; NEON/SVE-128: 8); no gathers |
| blend on t ≥ Shifter + 4 | j ≥ 16 on the integer index (equivalent for every lane whose table value is used) |
| `vminpd {sae}`, `vcmppd`, masked FMAs, `vblendmpd` | `Min`, `Lt`, `IfThenElse` |
| NaN propagation | NaN lanes are set to `x + x` (the quieted input, as the SVML code returns) |
| `vrcp14pd` | the instruction on AVX-512 targets; elsewhere its bit-exact emulation (below) |

The result depends on the `vrcp14pd` bits (Rcp₂ keeps the rounding error of RN(Rcp(1 + eps)), so R and
the final rounding change with Rcp). The instruction is defined bit-exactly by Intel's reference
RECIP14.c, which `common/include/hwymath/hwy/avx512_emul-inl.h` emulates without gathers (see asin's README; its
domain, positive normal inputs with a normal reciprocal, covers Y ∈ [1, 2¹²⁸]):

- the results are bit-identical to Intel's AVX-512 code (`svml_atan8_ha`) on every FMA target (AVX3_SPR,
  AVX3_DL, AVX3, AVX2 under SDE; AArch64 NEON/SVE/SVE2 identical to x86-64 under QEMU);
- the emulation costs time: AVX2 6.3 ns/value on the Core i5-9500 (4.5 ns with the earlier,
  non-identical trunc16(1/Y) substitute).

Targets without FMA (SSE2/SSSE3/SSE4) use `std::atan`. Intel's libsvml AVX2 code (`__svml_atan4_ha_l9`)
is a different, less accurate algorithm (0.73 ulp); the Highway AVX2 version is ~1.5x slower and keeps the
AVX-512 algorithm's 0.57 ulp.

## Tests (`make atan-test`, `make atan-avx512-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_datan_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,004,201 inputs, 0 mismatches |
| `test_atan8` | Intel's `__svml_atan8_ha`, all 8 lanes bit-exact | 2,004,208 inputs, 0 mismatches (Intel SDE `-spr`); skipped (77) without AVX-512 |
| `hwy_test` | every Highway target: error bound (max ≤ 0.6 ulp, none > 1 ulp), special inputs (NaN, ±inf, ±0, denormals), identical results on all FMA targets; with AVX-512: bit-identity of every FMA target with the C translation `svml_atan8_ha` | all checks pass (AVX3_SPR, AVX3_DL, AVX3, AVX2 under SDE; AVX2 natively) |
| `make atan-insns` | arithmetic instruction mix of the compiled vector path vs the assembly | identical |
| `make atan-xcheck` (`pixi run -e cross`) | Highway NEON, SVE and SVE2 under QEMU (SVE at 128, 256, 512, 1024 and 2048 bits, i.e. 2 to 32 lanes) vs x86-64 AVX2: inputs and results bit for bit (`common/tools/xcheck.sh` with `functions/asin/tools/xcheck_invtrig.cc`) | identical on all 11 AArch64 configurations (200,000 random inputs plus the edge cases) |

Inputs (`tests/inputs.h`): special values (±0, ±inf, NaNs with payloads, denormals, ±DBL_MAX, 2¹²⁸, ...),
±50 ulps around every threshold of the vector code (7.875, 3.875…4.125, 2¹²⁸) and of the scalar routine
(0.125, 2⁵⁴, 2⁻⁵³, 2⁻¹⁰²², ...), ±3 ulps around every table point j/4 and interval end (2j+1)/8, random
values (uniform in [−8, 8], log-uniform magnitudes in [2⁻⁶⁰, 2⁷⁰], arbitrary bit patterns).
`make atan-hashes` writes the result hashes (`svml` per target, `svml-c svml`) for the cross-machine
check; the exhaustive check of the `vrcp14pd` emulation (`make common-test`) runs in `make test`
and `atan-avx512-test`.
