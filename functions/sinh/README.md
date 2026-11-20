# sinh (double precision)

Source studied: [`svml_z0_sinh_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_sinh_d_ha.s)
(AVX-512, high accuracy; entry `__svml_sinh8_ha`, callout `__svml_dsinh_ha_cout_rare_internal`). Accuracy
and speed: [Performance.md](Performance.md). cosh and tanh are in [`../cosh`](../cosh) and [`../tanh`](../tanh).

| path | contents |
|---|---|
| `svml/svml_sinh_d_ha.c` | bit-identical C/AVX-512 translation of the vector routine and of the scalar callout |
| `svml/svml_sinh_d_ha_constants.inc` | constants and tables, generated bit-exactly from the assembly by `tools/gen_svml_constants.py` |
| `hwy/` | Google Highway version of the same algorithm, bit-identical on every target with FMA |
| `tests/` | `test_rare.c`, `test_sinh8.c`, `hwy_test.cc`, `hwy_xcheck.cc`, test inputs (`inputs.h`) |
| `tools/` | `gen_svml_constants.py`, `xcheck.sh` (AArch64 NEON/SVE under QEMU) |
| `sinh_d.h`, `function.mk` | declarations; build rules and the `sinh-*` make targets |

## Algorithm (as read from the assembly)

Vector path, 8 lanes, a = |x|, sign or-ed back at the end:

1. m = a·log2(e) + 1.5·2^48 (FMA, round to nearest): N = m − 1.5·2^48 = M + j/16 is a·log2(e) rounded to a
   multiple of 1/16, j = low 4 bits of m. If m is exactly 1.5·2^48 + 1/16 (N = 1/16, i.e. a in about
   [ln2/32, 3ln2/32)), N = 0 is used instead (`vcmppd` + `vblendmpd`), so small arguments use r = a directly.
2. r = a − N·ln2_hi − N·ln2_lo (two FMAs), |r| ≤ ln2/32.
3. Four 16-entry tables (`vpermt2pd`): Tp = 2^(j/16)/2 = TpH + TpL and Tm = 2^(−j/16)/2 = TmH + TmL. They
   are scaled by 2^±M by adding/subtracting M<<52 to/from their bit patterns (integer `vpaddq`/`vpsubq`);
   TmL is set to 0 where |TmL| < M<<52 (compared as doubles), i.e. where the subtraction would underflow.
4. sinh(a) = (Tp − Tm)·cosh(r) + (Tp + Tm)·sinh(r) with D = TpH − TmH, S = TpH + TmH:
   lo = (TpL − (rounding error of D)) − TmL; result = D + [D·r²(c2 + c4r² + c6r⁴ + c8r⁶) + lo +
   S·(r + r³(s3 + s5r² + s7r⁴))] (Horner FMAs, two independent chains for the even and odd parts).
5. Lanes whose high word of |x| exceeds 0x4084ee33 (|x| ≥ 669.78, inf, NaN) are recomputed by the scalar
   callout, one lane at a time.

Scalar callout (64-entry table of 2^(j/64) as {hi, lo} pairs, plus the entry 2.0 used for 2^(−0/64)·2):

| \|x\| | result |
|---|---|
| NaN, ±inf | x + x |
| 0, denormal | x·(2^−1022(1 + 2^−52)) + x (raises underflow/inexact for denormals) |
| < 2^−53 | (1 + tiny)·x |
| < 0.0272 | x + x³(s3 + s5x² + s7x⁴ + s9x⁶), Horner |
| < 21.49 | (e^a − e^−a)/2 in double-double arithmetic: Tp = 2^(M−1)·2^(j/64), Tm = 2^(−M−1)·2^(−j/64), r = rh + rl, Dekker splits (2^27 + 1) of rh and of Tp + Tm, ~40 operations |
| < 710.476 | e^a/2 = 2^(M−1)·2^(j/64)·(1 + p(r)), degree-6 polynomial |
| ≥ 710.476 | DBL_MAX·x = ±inf, raises overflow, return code 3 |

The sign is applied by overwriting the top bit of the result. The vector path calls the callout only for
|x| ≥ 669.78, i.e. the last two rows and NaN/inf.

## The C translation (`svml/`)

`svml_sinh8_ha(__m512d)` keeps every instruction, constant and rounding mode (`{rn-sae}`/`{sae}`) of the
assembly; `svml_sinh_ha_array(x, y, n)` loops with a masked tail. Two `asm("" : "+v"(v))` optimization
barriers stop GCC from merging the sign/index bitwise operations into `vandnpd`/`vpternlogq`, so that the
arithmetic instruction mix matches the assembly exactly (`make sinh-insns`); they do not affect results.
`svml_dsinh_ha_cout_rare()` is plain C (no FMA contraction), statement by statement in the order of the
assembly; values that only exist to raise flags are `volatile`.

## Highway version (`hwy/hwy_sinh.cc`)

Run-time dispatch over all Highway targets; `hwy_sinh_array(x, y, n)`.

- **Targets with FMA (`HWY_NATIVE_FMA`): the SVML algorithm, bit-identical to Intel's AVX-512 code.** The
  vector path only uses round-to-nearest arithmetic, so nothing needs a rounding-mode emulation. The
  AVX-512-specific operations are replaced by exact equivalents:

  | AVX-512 | Highway |
  |---|---|
  | `vpermt2pd` (16-entry tables) | see the table lookups below |
  | `vcmppd {sae}` into `k` + `vblendmpd` / masked `vxorpd` | `Eq`/`Lt` on doubles (same ordered, quiet predicates) + `IfThenElse` / `IfThenZeroElse` |
  | `vpsrlq` + `vpmovqd` + `vpcmpgtd` + `vmovmskps` (high word of \|x\| > 0x4084ee33) | one 64-bit signed compare: bits(\|x\|) > 0x4084ee33ffffffff |
  | integer exponent arithmetic | `ShiftLeft`, `Add`/`Sub` on `BitCast` u64 vectors |

  Callout lanes are recomputed with the C callout (`svml_dsinh_ha_cout_rare`) in an out-of-line, `cold`
  loop.
- Targets without FMA (SSE2/SSSE3/SSE4): `std::sinh` (the FMA results cannot be reproduced exactly).

### Table lookups without gathers

The four tables are looked up with the same index. Per vector length:

| lanes | lookup |
|---|---|
| 8 (AVX-512, SVE-512) | one `TwoTablesLookupLanes` per table = Intel's `vpermt2pd` |
| ≥ 16 (SVE/RVV ≥ 1024 bit) | `TableLookupLanes` of the 16 entries (`LoadN`) |
| 4 on x86 (AVX2) | **by rows**: the tables are interleaved at compile time into rows {TpH, TpL, TmH, TmL}[j]; the index vector is stored, and each lane's row is loaded as two 128-bit halves that are combined (`Combine`, i.e. `vinsertf128` from memory) and interleaved: 8 half-row loads + 4 inserts + 4 unpacks for all four tables |
| 2 (NEON, SVE-128) | by rows: two 128-bit loads per lane, 4 `InterleaveLower/Upper` for all four tables |
| 4 elsewhere (SVE-256) | permute/blend tree: per table 2 `TwoTablesLookupLanes` + 1 `IfThenElse` on bit 3 |

Measured on the Core i5-9500 (AVX2 target, ns/value, several alternating runs on a shared machine), with
`-DHWYSINH_ROWS=0/1/2`:

| lookup on AVX2 | sinh | cosh | tanh |
|---|---|---|---|
| permute/blend tree (`TwoTablesLookupLanes` ×2 + `IfThenElse` per table) | 2.41–2.45 | 1.72–1.74 | 5.92–5.95 |
| rows, full-row loads + 4×4 transpose (4 loads + 8 shuffles per 4 tables) | 2.02–2.12 | 1.33–1.37 | 3.34–3.50 |
| **rows, half-row loads + inserts + unpacks (default)** | 1.97–2.07 | 1.34–1.38 | 3.25–3.36 |

The tree costs 7 shuffle/blend instructions per table on AVX2 (Highway's `TwoTablesLookupLanes` is two
`vpermps` and a `vblendvps`), all on port 5; the row loads move the work to the load ports. Gathers were not
considered (very slow with the GDS microcode mitigation). The 2- and 4-lane paths are also exercised on
wider targets by `hwy_test` (vectors capped at 2 and 4 lanes) and on AArch64 by `make sinh-xcheck`.

Deviations from Intel's algorithm: none (bit-identical everywhere it runs with FMA).

## Tests

| test | checks | result (this machine, i5-9500; AVX-512 under Intel SDE `-spr`) |
|---|---|---|
| `test_rare` (`make sinh-test`) | C callout vs Intel's `__svml_dsinh_ha_cout_rare_internal` (made global with `sed`): results, return codes, FP exception flags | 20,007,045 inputs, 0 mismatches |
| `test_sinh8` (`make sinh-avx512-test`) | `svml_sinh8_ha` vs `__svml_sinh8_ha`, all lanes, bit for bit | 2,007,048 inputs, 0 mismatches |
| `make sinh-insns` | arithmetic instruction mix (mnemonic + rounding) of the compiled vector path vs the assembly | identical |
| `hwy_test` (`make sinh-test`, `sinh-avx512-test`) | every compiled+supported Highway target: max error < 1 ulp vs `sinhl`, special values exact, all SVML-algorithm targets identical (also with vectors capped at 2/4 lanes); with AVX-512: bit-identical to `svml_sinh_ha_array` | AVX3_SPR, AVX3_DL, AVX3, AVX2: bit-identical to the C translation (2,007,045 inputs, SDE) |
| `make sinh-hashes` | `arch sinh target svml n hash` lines; the C translation's hash appears as target `svml-c` when run with AVX-512 | AVX2 hash = `svml-c` hash |
| `make sinh-xcheck` (pixi `cross` env) | static NEON, SVE and SVE2 builds under QEMU (SVE at 128/256/512/2048 bits) vs x86-64 AVX2, bit for bit | all bit-identical (300,000 random inputs) |

Inputs (`tests/inputs.h`, deterministic): ±0, ±inf, NaNs, denormals, DBL_MAX; ±50 doubles around every
threshold (2^−53, 0.0272, 21.49, 710.476 = log(2·DBL_MAX), the vector callout 669.78, the N = 1/16 range
ln2/32 and 3ln2/32, ...), table boundaries (k + 1/2)·ln2/64, uniform values in [−720, 720] and [−25, 25],
±m·2^e with e in [−60, 3] (the sinh(x) ≈ x region), and random bit patterns.
