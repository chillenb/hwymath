# tan (double precision)

Source studied: [`svml_z0_tan_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_tan_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md).

tan shares the large-argument reduction with sin and cos; the shared parts (reduction in C and Highway,
table lookups, test templates and inputs, `trig_common.mk`) live in [`functions/sin`](../sin/README.md).

| path | contents |
|---|---|
| `svml/svml_tan_d_ha.c` | bit-identical C/AVX-512 translation of the assembly and its scalar callout |
| `svml/svml_tan_d_ha_constants.inc` | its constants, generated from the assembly by `tools/gen_svml_constants.py` (which also checks that the reduction data equal sin's) |
| `hwy/hwy_tan.cc` | Highway version |
| `hwy/rcp14_hwy-inl.h` | exact Highway emulation of `vrcp14pd` |
| `tests/` | `test_rare.c`, `test_tan8.c`, `test_rcp14.c`, `rcp14_inputs.h`, `hwy_test.cc`, `extras.c` |
| `tan_d.h` | declarations of the C translation |

## Algorithm

**Main-path reduction** (|x| <= 65536.0035): N = round(x*16/pi) (`y = x*16/pi + 1.5*2^52`), j = N mod 16,
R + Rlo = x - N*(PI1 + PI2 + PI3) (|R| <= pi/32).

**Large-path reduction** (|x| > 65536.0035, +-inf, NaN; all 8 lanes recomputed if one needs it): the
Payne-Hanek reduction of sin/cos (`functions/sin/svml/svml_trig_reduction.h`), applied to x with its
sign, gives x = M*pi/256 + r + rlo (M mod 512). Then M = 16k + i with k = (M + 8 - sign(r)) >> 4 (round
to nearest, ties by the sign of r; 32-bit integer operations and `vcvtdq2pd` in the assembly), |i| <= 8,
j = k mod 16 and R + Rlo = r + i*pi/256 + rlo (pi/256 in two parts). The large path jumps back into the
evaluation of the main path.

**Evaluation**: with Th + Tl = tan(j*pi/16) from two 16-entry tables (`vpermt2pd`; for j = 8, the pole,
Th = -DBL_MAX and Tl = -2^971, which turns the quotient below into -1/tan R):
- P + Plo = tan R = R + R^3*poly(R^2) + Rlo (degree 4 in R^2), split exactly;
- numerator Num + Nlo = Th + P [+ Tl + Plo] (two-sum);
- denominator D - Dlo = 1 - Th*P [- Th*Plo - Tl*P], with the exact error of the FMA;
- tan x = (Num + Nlo)/(D - Dlo) without a division: rcp = `vrcp14pd`(D) (14-bit approximation), one
  Newton-Raphson step on D - Dlo, Q = Num*rcp, and one correction of the quotient,
  result = Q - rcp*((Q*D - Num) - Q*Dlo - Nlo).

**Callout**: tan(+-inf) = NaN (invalid, return code 1), tan(NaN) = x*x; same routine as sin's.

## vrcp14pd: exact emulation

The result depends on the exact bits of `vrcp14pd`: replacing it by a correctly rounded reciprocal
changes 0.20 % of the results (0.12 % on the main path), and changing its result by one unit of its 17th
significant bit changes 0.50 % (measured under SDE on the 2M test inputs). Exact emulation is possible
because Intel documents the instruction's results: the SDM refers to Intel's "Reference Implementations
for IA Approximation Instructions VRCP14, VRSQRT14, VRCP28, VRSQRT28, and VEXP2" (`recip14.c`, BSD-3),
which is stated to be numerically exact. For a significand 1.f (not a power of 2) it reads only the top
16 fraction bits: with i = the top 6 and l = the next 10,
rcp14 = truncate_to_17_bits(a[i] - b[i]*(m17 - y[i])) / 2^18, m17 = 1.f truncated to 16 fraction bits,
y[i] = 1 + 1/128 + i/64, and a 64-entry table of 18-bit a[i] and 10-bit b[i].

**Why the emulation is exact.** Substituting m17 - y[i] = (l - 512)/2^16, the value before truncation is
(256 a[i] + 512 b[i] - b[i] l)/(256 * 2^18) with an integer numerator; over all 64 x 1024 (i, l) pairs
that numerator divided by 512 lies in [65536.006, 131068.5] = [2^16, 2^17), so the exponent of the
intermediate value never changes and truncation to 17 significant bits is `floor(numerator / 512)`:
rcp14 = floor((256 a[i] + 512 b[i] - b[i] l) / 512) * 2^-17, an integer computation. The remaining cases
of the reference are exponent manipulations: biased exponent 2045 - e (2046 - e for powers of 2, which
are exact), results with e = 2045/2046 shifted into the denormal range (by truncation), denormal inputs
normalized first (|x| <= 2^-1024 gives +-inf), +-0 -> +-inf, +-inf -> +-0, NaN -> quiet NaN (DAZ = FTZ = 0,
SVML's environment). `common/include/hwymath/rcp14_ref.h` implements exactly this; `hwy/rcp14_hwy-inl.h` implements the same
integer formula in vectors: the pairs {a[i], b[i]} are packed in 28-bit words and looked up with permutes
(`Lookup<64>` on 32-bit lanes: 2 `vpermt2d` + 1 select on AVX-512, 4 x (2 `vpermd` + blend) + 3 selects on
AVX2; no gathers), b[i]*l with `MulEven`, denormal inputs by scaling with 2^52, the special cases by
selects.

**Verification.** `test_rcp14` compares the reference with `vrcp14pd` (here: Intel SDE `-spr`) and
`hwy_test` compares the Highway emulation with the reference on every target (AVX-512 targets use
`vrcp14pd` itself, so there it checks the instruction), on 2,474,579 inputs
(`tests/rcp14_inputs.h`): all 2^16 leading fraction patterns at 10 exponents including those with special
handling (1, 2, 0x3fe, 0x3ff, 0x400, 0x7fc, 0x7fd, 0x7fe), both signs, with and without random ignored
bits; denormals with every leading pattern after normalization; every exponent with random fractions;
special values; 1M random bit patterns: 0 mismatches each. Since the instruction reads only the sign, the
exponent, 16 fraction bits and whether the remaining fraction is zero (power of 2), and the exponent
enters only through the mirrored-exponent formula and the boundary cases tested above, this covers its
whole behaviour. Not checked here on real AVX-512 hardware (CI runs `test_rcp14` and `test_tan8` on a
Xeon); whether AMD's Zen 4/5 `vrcp14pd` matches Intel's reference is unknown: if not, Intel's own
`__svml_tan8_ha` gives different results there, while the Highway version keeps reproducing the Intel
results.

**Main path: a cheaper variant.** On the main path the denominator D = 1 - Th*P satisfies
|P| <= tan(pi/32 + tiny) < 0.0986: for j != 8, |Th| <= 5.03, so D is in [0.5, 1.5]; for j = 8,
|D| <= 1 + DBL_MAX*0.0986 < 2^1021, and |D| >= DBL_MAX*|P| - 1 > 2^960 because no double comes closer
than about 2^-62 to a multiple of pi/2 (the closest one, 6381956970095103*2^797, is at 2^-60.9). So D is always normal with biased exponent in [0x3fe, 0x7fb], and inputs with
x = inf/NaN (D = NaN) get the callout. The main loop therefore uses `Rcp14Normal`, the same formula
without the special cases (exponents 1..0x7fc; powers of 2 included: D = 1 exactly for j = 0), and the
rare out-of-line path uses the complete emulation. Both are tested against the reference (the normal
variant on the 2,079,365 inputs with exponents in its range).

## The C translation (`svml/svml_tan_d_ha.c`)

`svml_tan8_ha(__m512d)`, `svml_tan_ha_array()`, `svml_dtan_ha_cout_rare()`, as for sin (see
[its README](../sin/README.md) for the differences from the assembly that do not affect results). The
evaluation is a separate `noinline` function `svml_tan8_ha_eval` because the assembly reaches it from both
reduction paths (`.LBL_1_2`); `make tan-insns` compares the reduction part and the evaluation part
separately (the build inserts a global label `__svml_tan8_ha_eval` at `.LBL_1_2` in the copy of the
assembly that is assembled). Intel's libsvml (`__svml_tan8_ha_z0`) agrees with the translation on all
test inputs except tan(+-inf) (+NaN instead of -NaN).

## Highway version (`hwy/hwy_tan.cc`)

The SVML algorithm on every target with FMA, bit-identical to Intel's AVX-512 code; `std::tan` on targets
without FMA. Emulations: as for sin (`{rn-sae}` arithmetic, `vpermt2pd` of the two 16-entry tables,
masks, `vpmullq`, the gathers of the large path, `vcvtdq2pd`/32-bit integer operations of the large
path by 64-bit integer operations and `ConvertTo`, all exact) plus `vrcp14pd` as above. **Deviations:
none** on FMA targets. Bit-identity verified under SDE for AVX3_SPR, AVX3_DL, AVX3 and AVX2 against the
C translation, natively via the AVX2 hash, and by `make tan-xcheck` for NEON_WITHOUT_AES, NEON, SVE and
SVE2 (QEMU, SVE at 128 to 2048 bits) against x86-64. The `vfnmsub` of P = tan R uses
`NegMulAdd(a, b, Neg(c))`: Highway's `NegMulSub` is `Neg(MulAdd(a, b, c))` on NEON, which made
tan(-0) = +0 there (found by CI on a Neoverse-N2).

Cost: the exact `vrcp14pd` emulation is the expensive part on AVX2 (about 2.4 of the ~5 ns/value; the
64-entry lookup needs 8 `vpermd` per vector of 4 doubles). With a hardware division instead, the
algorithm would take ~2.6 ns/value (and be slightly more accurate: 0.513 instead of 0.568 ulp), but would
no longer reproduce Intel's results.

## Tests (`make tan-test`, `make tan-avx512-test RUN="$SDE -spr --"`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's `__svml_dtan_ha_cout_rare_internal`: results, whether written, return codes, FP flags | 20,035,562 inputs, 0 mismatches |
| `test_rcp14` | `vrcp14pd` vs `common/include/hwymath/rcp14_ref.h` (Intel SDE here) | 2,474,584 inputs, 0 mismatches |
| `test_tan8` | Intel's `__svml_tan8_ha`, bit for bit (Intel SDE) | 20,035,568 inputs (10,783,279 in the large path), 0 mismatches |
| `hwy_test` | every Highway target: error < 1 ulp, special values, agreement between targets, bit-identity with the C translation (SDE: AVX3_SPR, AVX3_DL, AVX3, AVX2), `vrcp14pd` emulations vs the reference | all checks pass |
| `make tan-insns` | instruction mix of the reduction and of the evaluation vs the assembly | identical |

Inputs: `functions/sin/tests/inputs.h` (described in sin's README; it includes tan's threshold
65536.0035 and the doubles next to multiples of pi/2, i.e. the poles and zeros of tan).
