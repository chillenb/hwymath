# cosh (double precision)

Source studied: [`svml_z0_cosh_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_cosh_d_ha.s)
(AVX-512, high accuracy; entry `__svml_cosh8_ha`, callout `__svml_dcosh_ha_cout_rare_internal`). Accuracy
and speed: [Performance.md](Performance.md). The structure follows [`../sinh`](../sinh); see its README for
the table-lookup measurements, which cover all three hyperbolic functions.

| path | contents |
|---|---|
| `svml/svml_cosh_d_ha.c` | bit-identical C/AVX-512 translation of the vector routine and of the scalar callout |
| `svml/svml_cosh_d_ha_constants.inc` | constants and tables, generated bit-exactly from the assembly by `tools/gen_svml_constants.py` |
| `hwy/` | Google Highway version of the same algorithm, bit-identical on every target with FMA |
| `tests/` | `test_rare.c`, `test_cosh8.c`, `hwy_test.cc`, `hwy_xcheck.cc`, test inputs (`inputs.h`) |
| `tools/` | `gen_svml_constants.py`, `xcheck.sh` (AArch64 NEON/SVE under QEMU) |
| `cosh_d.h`, `function.mk` | declarations; build rules and the `cosh-*` make targets |

Why a separate directory from sinh: the two share the idea (2^(±j/16) tables, exponent arithmetic) and
the 16-entry lookup helper, but little code. sinh has the N = 1/16 → 0 adjustment, a fourth table (TmL)
with an underflow guard, sign handling, a different reconstruction, and a scalar callout with a Taylor
branch and a ~40-operation double-double path; cosh has none of these. Merging them would have meant
two kernels and two callouts in one file. The common part (the lookups, ~60 lines) is duplicated; a shared
helper is proposed in the porting report.

## Algorithm (as read from the assembly)

Vector path, 8 lanes, a = |x|:

1. m = a·log2(e) + 1.5·2^48 (FMA, round to nearest): N = M + j/16, j = low 4 bits of m (no N = 1/16
   adjustment, unlike sinh).
2. r = a − N·ln2_hi − N·ln2_lo, |r| ≤ ln2/32.
3. Three 16-entry tables (`vpermt2pd`, indexed directly by m): Tp = 2^(j/16)/2 = TpH + TpL and
   TmH ≈ 2^(−j/16)/2, scaled by 2^±M by integer addition/subtraction of (m<<48) & 0x7ff0000000000000.
4. cosh(a) = (Tp + Tm)·cosh(r) + (Tp − Tm)·sinh(r):
   result = TpH + (TmH + ((Tp − Tm)·(r + r³(s3 + s5r² + s7r⁴)) + ((Tp + Tm)·r²(c2 + c4r² + c6r⁴) + TpL))).
5. Lanes whose high word of |x| exceeds 0x40861d99 (|x| ≥ 707.70, inf, NaN) are recomputed by the scalar
   callout.

Scalar callout (64-entry table of 2^(j/64), plus 2.0 for j = 64):

| \|x\| | result |
|---|---|
| NaN, ±inf | x·x |
| < 2^−53 (incl. 0, denormals) | 1 + \|x\| |
| < 21.49 | Tp = 2^(M−1)·2^(j/64), Tm = 2^(−M−1)·2^(−j/64): S + ((TpL + (TpL − TmL)·r + TmL + D·(sinh(r) − r) + r·D) + S·(cosh(r) − 1) + low part of S), D = TpH − TmH, S = TpH + TmH |
| < 710.476 | e^a/2, degree-6 polynomial |
| ≥ 710.476 | DBL_MAX·DBL_MAX = +inf, raises overflow, return code 3 |

## The C translation (`svml/`)

`svml_cosh8_ha(__m512d)` keeps every instruction, constant and rounding mode of the assembly;
`svml_cosh_ha_array(x, y, n)` loops with a masked tail; `svml_dcosh_ha_cout_rare()` is plain C in the
assembly's order (no FMA contraction; the overflow product is `volatile`). Instruction mix identical to
the assembly without any optimization barrier (`make cosh-insns`).

## Highway version (`hwy/hwy_cosh.cc`)

Same design as sinh: on targets with FMA the SVML algorithm, bit-identical to Intel's AVX-512 code (only
round-to-nearest arithmetic; `vpermt2pd` → lookups below; `vpmovqd`/`vpcmpgtd`/`vmovmskps` on the high
words → one 64-bit compare bits(|x|) > 0x40861d99ffffffff; exponent arithmetic on `BitCast` u64 vectors);
callout lanes use the C callout; `std::cosh` on targets without FMA.

Lookups of the three tables: `TwoTablesLookupLanes` (= `vpermt2pd`) on 8 lanes, `TableLookupLanes` on ≥ 16,
rows {TpH, TpL, TmH, 0}[j] on AVX2 (6 half-row loads + 3 inserts + 3 unpacks) and on 2-lane targets (4 loads
+ 3 zips), and a permute/blend tree on other 4-lane targets (SVE-256). On the i5-9500 the row lookup runs at
1.33–1.38 ns/value against 1.72–1.74 for the tree (details: [`../sinh/README.md`](../sinh/README.md)).

Deviations from Intel's algorithm: none.

## Tests

| test | checks | result (this machine, i5-9500; AVX-512 under Intel SDE `-spr`) |
|---|---|---|
| `test_rare` (`make cosh-test`) | C callout vs `__svml_dcosh_ha_cout_rare_internal`: results, return codes, FP flags | 20,006,641 inputs, 0 mismatches |
| `test_cosh8` (`make cosh-avx512-test`) | `svml_cosh8_ha` vs `__svml_cosh8_ha`, bit for bit | 2,006,648 inputs, 0 mismatches |
| `make cosh-insns` | arithmetic instruction mix vs the assembly | identical |
| `hwy_test` | every Highway target: < 1 ulp vs `coshl`, special values, identical targets (also capped at 2/4 lanes), bit-identical to the C translation with AVX-512 | AVX3_SPR, AVX3_DL, AVX3, AVX2: bit-identical to the C translation (2,006,641 inputs, SDE) |
| `make cosh-hashes` | hash lines, `svml-c` = C translation | AVX2 hash = `svml-c` hash |
| `make cosh-xcheck` | NEON, SVE/SVE2 at 128–2048 bits under QEMU vs x86-64 | all bit-identical (300,000 random inputs) |

Inputs (`tests/inputs.h`): special values, ±50 doubles around 2^−53, 21.49, 710.476, the vector callout
threshold 707.70, ln2/64, ln2/32, ...; table boundaries (k + 1/2)·ln2/64; uniform values in [−720, 720]
and [−25, 25]; ±m·2^e with e in [−60, 3]; random bit patterns.
