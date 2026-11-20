# tanh (double precision)

Source studied: [`svml_z0_tanh_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_tanh_d_ha.s)
(AVX-512, high accuracy; entry `__svml_tanh8_ha`, callout `__svml_dtanh_ha_cout_rare_internal`). Accuracy
and speed: [Performance.md](Performance.md). Same layout as [`../sinh`](../sinh) and [`../cosh`](../cosh)
(the algorithm is unrelated to theirs: piecewise polynomials, no exponential).

| path | contents |
|---|---|
| `svml/svml_tanh_d_ha.c` | bit-identical C/AVX-512 translation of the vector routine and of the scalar callout |
| `svml/svml_tanh_d_ha_constants.inc` | 19 coefficient tables and masks, generated bit-exactly from the assembly by `tools/gen_svml_constants.py` |
| `hwy/` | Google Highway version of the same algorithm, bit-identical on every target with FMA |
| `tests/` | `test_rare.c`, `test_tanh8.c`, `hwy_test.cc`, `hwy_xcheck.cc`, test inputs (`inputs.h`) |
| `tools/` | `gen_svml_constants.py`, `xcheck.sh` (AArch64 NEON/SVE under QEMU) |
| `tanh_d.h`, `function.mk` | declarations; build rules and the `tanh-*` make targets |

## Algorithm (as read from the assembly)

Vector path, 8 lanes:

1. Subinterval index from the high words (32-bit integer operations on `vpmovqd`-narrowed lanes):
   j = clamp(((hi & 0x7ff80000) − 0x3fc00000) >> 19, 0, 15) — j = 0 for |x| < 0.1875, then half-binades
   [0.1875, 0.25), [0.25, 0.375), ..., [16, 24), and j = 15 for |x| ≥ 24.
2. y = |x| − B[j], B = midpoint of the subinterval (0 for j = 0 and j = 15).
3. 19 tables of 16 entries, all looked up with `vpermt2pd`: B, P0 (high and low part), P1 ... P16.
   tanh|x| = P0h + (P1·y + (y·(y·(P2 + P3y + ... + P16y¹⁴)) + P0l)): a 15-step Horner FMA chain. For
   j = 15 the tables hold P0 = 1 and zeros (saturation: tanh = 1 for |x| ≥ 24, rounding to 1 already
   happens from about 19.06). For j = 0, B = 0 and P1 = 1: tanh(x) = x + O(x³), exact for tiny x.
4. Result = tanh|x| with the sign bit of x or-ed in (tanh(−0) = −0).
5. Lanes with (hi & 0x7ff80000) > 0x7fe00000, i.e. |x| ≥ 1.5·2^1023, inf and NaN, go to the callout.

Scalar callout: NaN → x + x; ±inf and biased exponent 2046 → ±1 (from a two-entry table); any other
input → x·(1 + x) (never reached from the vector path; translated as is).

## The C translation (`svml/`)

`svml_tanh8_ha(__m512d)` keeps every instruction, constant and rounding mode of the assembly, including the
32-bit index arithmetic on `__m256i` (`vpand`, `vpsubd`, `vpmaxsd`, `vpminsd`, `vpsrld`, `vpmovzxdq`);
`svml_tanh_ha_array(x, y, n)` loops with a masked tail; `svml_dtanh_ha_cout_rare()` is plain C. One
`asm("" : "+v"(v))` barrier keeps GCC from folding the sign `vandpd`/`vorpd` into a `vpternlogq`, and the
zero for `vpmaxsd` is produced with the assembly's `vxorps` (GCC would emit `vpxor`), so that the
arithmetic instruction mix matches (`make tanh-insns`); neither affects results.

## Highway version (`hwy/hwy_tanh.cc`)

On targets with FMA the SVML algorithm, bit-identical to Intel's AVX-512 code; `std::tanh` on targets
without FMA. Emulation of the AVX-512-specific operations (all exact; the arithmetic is round-to-nearest
only):

| AVX-512 | Highway |
|---|---|
| `vpermt2pd` ×19 | see below |
| `vpsrlq 32` + `vpmovqd` + 32-bit and/sub/max/min/shift + `vpmovzxdq` (index) | u = (bits(x) >> 51) & 0xfff in 64-bit lanes, then `Max`/`Min`/`Sub` with 0x7f8/0x807 in 32-bit lanes (`Repartition<int32_t>`): j = clamp(u − 0x7f8, 0, 15); the zero upper halves of the 64-bit lanes stay zero (max(0, 0x7f8) − 0x7f8) |
| `vpcmpgtd` + `vmovmskps` (callout lanes) | 64-bit signed compare bits(\|x\|) > 0x7fe7ffffffffffff |

Table lookups (19 tables, one index), per vector length:

| lanes | lookup |
|---|---|
| 8 (AVX-512, SVE-512) | `TwoTablesLookupLanes` per table = Intel's `vpermt2pd` |
| ≥ 16 | `TableLookupLanes` |
| 4 on x86 (AVX2) | by rows: the tables are interleaved at compile time into rows of 20 doubles {B, P0h, P0l, P1, ..., P16, 0}[j]; for each group of 4 tables: 8 half-row loads, 4 `vinsertf128`, 4 unpacks (5 groups) |
| 2 (NEON, SVE-128) | by rows: per group of 4 tables 4 128-bit loads + 4 zips |
| 4 elsewhere (SVE-256) | permute/blend tree (2 `TwoTablesLookupLanes` + 1 `IfThenElse` per table) |

The groups are fetched in the order the Horner chain consumes them, so that only the current group is
live. On the i5-9500 (AVX2), the row lookups run at 3.25–3.36 ns/value (3.34–3.50 with full-row loads and a
4×4 transpose) against 5.92–5.95 for the permute/blend tree: 19 lookups × 7 port-5 instructions dominate
the tree version, while the row version is bound by the Horner chain (16 dependent FMAs). Details of the
measurement: [`../sinh/README.md`](../sinh/README.md).

Deviations from Intel's algorithm: none.

## Tests

| test | checks | result (this machine, i5-9500; AVX-512 under Intel SDE `-spr`) |
|---|---|---|
| `test_rare` (`make tanh-test`) | C callout vs `__svml_dtanh_ha_cout_rare_internal`: results, return codes, FP flags | 20,005,683 inputs, 0 mismatches |
| `test_tanh8` (`make tanh-avx512-test`) | `svml_tanh8_ha` vs `__svml_tanh8_ha`, bit for bit | 2,005,688 inputs, 0 mismatches |
| `make tanh-insns` | arithmetic instruction mix vs the assembly | identical |
| `hwy_test` | every Highway target: < 1 ulp vs `tanhl`, special values (±0, ±inf → ±1, saturation, DBL_MAX), identical targets (also capped at 2/4 lanes), bit-identical to the C translation with AVX-512 | AVX3_SPR, AVX3_DL, AVX3, AVX2: bit-identical to the C translation (2,005,683 inputs, SDE) |
| `make tanh-hashes` | hash lines, `svml-c` = C translation | AVX2 hash = `svml-c` hash |
| `make tanh-xcheck` | NEON, SVE/SVE2 at 128–2048 bits under QEMU vs x86-64 | all bit-identical (300,000 random inputs) |

Inputs (`tests/inputs.h`): special values; ±50 doubles around every subinterval boundary 1.5·2^e and
2^(e+1) from 0.1875 to 24, the saturation points ~18.71 and ~19.06, the callout threshold 1.5·2^1023,
2^1023, 2^−26, 2^−54, 1e−300; uniform values in [−30, 30] and [−2, 2]; ±m·2^e with e in [−60, 5]; random
bit patterns.
