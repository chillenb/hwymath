# acos (double precision)

Source studied: [`svml_z0_acos_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_acos_d_ha.s)
(AVX-512, high accuracy; entry `__svml_acos8_ha`). Accuracy and speed: [Performance.md](Performance.md).
The algorithm shares its building blocks (and constants) with asin; see [../asin/README.md](../asin/README.md)
for the details common to both, in particular the mechanical callout translation and the `vrsqrt14pd`
deviation. acos has its own directory because the vector code (selection of the constants, signs, order
of the final additions), the callout and the tests differ; the shared parts are the generator
`functions/asin/tools/asm2c.py` and the machine model `functions/asin/svml/asm_model.h`.

| path | contents |
|---|---|
| `svml/svml_acos_d_ha.c` | bit-identical C/AVX-512 translation of the vector routine, and a `_portable` variant |
| `svml/svml_acos_d_ha_constants.inc` | vector-path constants, generated from the assembly by `tools/gen_svml_constants.py` |
| `svml/svml_dacos_ha_cout_rare.inc` | the scalar callout, translated mechanically by `functions/asin/tools/asm2c.py` |
| `hwy/` | portable Google Highway version |
| `tests/` | tests, benchmark, comparison implementations (Intel SVML, SLEEF), test inputs |
| `acos_d.h` | declarations of the C kernels |
| `function.mk` | build rules and the `acos-*` make targets |

## Algorithm

With a = |x|, Y = 0.5 − 0.5a (computed as 0.5 + 0.5·(−|x|)), R = min(a², Y), and P the degree-11 asin
polynomial (asin(s) = s + s³P(s²)):

- a < 0.5: acos(x) = π/2 − asin(x) = π/2_hi + (−(x + x·R·P(R)) + π/2_lo);
- a ≥ 0.5, x > 0: acos(x) = 2 asin(√Y) = s + s·R·P(R) with s = 2√Y = Sh − Sl;
- a ≥ 0.5, x < 0: acos(x) = π − 2 asin(√Y) = π_hi + (−(s + s·R·P(R)) + π_lo).

2√Y = Sh − Sl is the same double-double as in asin (from r = `vrsqrt14pd`(Y), set to 0 for Y < 2⁻²⁵⁵).
The pair (hi, lo) ∈ {(π/2_hi, π/2_lo), (π_hi, π_lo), (0, 0)} is selected with masked xor (zeroing) and
masked or; x < 0 is detected as "x ≤ R" (R ≥ 0). The sign of x is applied with xors before the final
addition of hi, so the low parts enter first. Lanes with |x| > 1 (−|x| < −1, including ±inf) are
recomputed by the scalar routine `__svml_dacos_ha_cout_rare_internal` (NaN, invalid exception, code 1; a
complete scalar acos otherwise). NaN inputs stay on the vector path, whose result is the quieted input
NaN **with its sign flipped** (the NaN passes through −|x| and the final sign xor).

## The C translation

`svml/svml_acos_d_ha.c` keeps the algorithm, the constants, the operation order and every instruction's
rounding mode: `svml_acos8_ha(__m512d)` (bit-identical to Intel's, `test_acos8`), `svml_acos_ha_array`,
`svml_acos8_ha_portable` (with the portable rsqrt approximation; the reference for the Highway targets)
and `svml_dacos_ha_cout_rare(&x, &r)` (the generated scalar callout, portable C). As in asin, the bitwise
operations use the AVX-512F integer forms and empty `asm` barriers keep GCC from fusing the sign
operations into `vpternlogq`, so that `make acos-insns` compares equal instruction mixes.

## Highway version (`hwy/`)

`hwy/hwy_acos.cc`: the same emulation as asin (IEEE round-to-nearest operations and FMAs for the
`{rn-sae}` instructions, `Min` for `vminpd`, `IfThenElse` for masked operations and blends; the masked
xor/or selection of the π constants becomes nested `IfThenElse`), NaN lanes set to −(x + x), lanes with
|x| > 1 through the C callout, and `vrsqrt14pd` replaced by **trunc16(RN(1/RN(√Y)))** on every FMA
target (see asin for why `vrsqrt14pd` cannot be reproduced and why the result depends on it):

- the results are bit-identical on every FMA target and machine, and equal to `svml_acos8_ha_portable`;
- they differ from Intel's in 483 of the 2,004,655 test inputs (0.024 %, |x| ≥ 0.5 only); max error 0.817
  ulp as Intel's, not correctly rounded 6.77 % vs 6.76 %;
- `hwy_acos_native_array` uses `vrsqrt14pd` on AVX-512 targets: bit-identical to Intel's code on the same
  CPU (verified under SDE), which proves that all other operations are emulated exactly.

Targets without FMA (SSE2/SSSE3/SSE4) use `std::acos`.

## Tests (`make acos-test`, `make acos-avx512-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dacos_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,004,655 inputs, 0 mismatches |
| `test_acos8` | Intel's `__svml_acos8_ha`, all 8 lanes bit-exact | 2,004,656 inputs, 0 mismatches (Intel SDE `-spr`); skipped (77) without AVX-512 |
| `hwy_test` | every Highway target: error bound (max ≤ 0.9 ulp, none > 1 ulp), special inputs (NaN, ±inf, \|x\| > 1, acos(1) = +0, acos(−1) = RN(π)), identical results on all FMA targets; with AVX-512: bit-identity of every FMA target with `svml_acos8_ha_portable` and of the native variant with `svml_acos8_ha` | all checks pass (AVX3_SPR, AVX3_DL, AVX3, AVX2 under SDE; AVX2 natively) |
| `make acos-insns` | arithmetic instruction mix of the compiled vector path vs the assembly | identical |

Inputs (`tests/inputs.h`, the same set as asin's): special values, ±50 ulps around every threshold of the
vector code and the callout, arguments approaching ±1, random values (uniform in [−1, 1] and
[−1.25, 1.25], log-uniform magnitudes, arbitrary bit patterns). `make acos-hashes` writes the result
hashes (`port`, `native`, `svml-c intel`, `svml-c port`) for the cross-machine check.
