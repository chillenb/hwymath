# atanh (double precision)

Source studied: [`svml_z0_atanh_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_atanh_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md). The two logarithms use the
same reduction and tables (with the opposite sign) as asinh's ([../asinh/README.md](../asinh/README.md)),
and the Highway version uses the AVX-512 emulations described there (`functions/asinh/hwy/avx512_emul-inl.h`).

| path | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_atanh_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py`) |
| `hwy/hwy_atanh.cc` | portable Google Highway version, bit-identical to Intel's code on every target with FMA |
| `tests/` | tests, comparison implementations (Intel SVML, SLEEF), test inputs |
| `atanh_d.h` | declarations of the C kernels |
| `function.mk` | build rules and the `atanh-*` make targets |

## Algorithm

With y = |x|, atanh(y) = (log(1 + y) − log(1 − y))/2, and the sign of x is ORed into the factor ½:

- Yp = 1 + y and Ym = 1 − y as double-doubles (Yp_low = y − (Yp − 1), Ym_low = y + (Ym − 1), so that
  Ym = Ym_high − Ym_low exactly);
- for each: Rcp = `vrcp14pd`(Y) rounded to 1+4 mantissa bits (integer add of 2⁴⁷, mask),
  R = Rcp·Y − 1 ± Rcp·Y_low, K = `vgetexppd`(Rcp), tables log(1 + k/16) (high/low, `vpermt2pd`/`vpermi2pd`),
  and a degree-8 polynomial; there is no separate small-input path (for tiny y, Rcp = 1 and Rp ≈ y, Rm ≈ −y);
- the difference is assembled as a double-double: high part (Km − Kp)·L2H + ThM − ThP + Rp − Rm, low part
  (Km − Kp)·L2L + TlM − TlP + Rp²·P(Rp) − Rm²·P(Rm) + the rounding errors of the high additions; then
  multiplied by ±½.
- Lanes with |x| ≥ 1 or NaN call the scalar routine: NaN → x·x (return code 0); ±1 → ±1/0 = ±inf with the
  divide-by-zero flag (code 2); |x| > 1 (and ±inf) → NaN with the invalid flag (code 1). The NaN is computed as
  0·inf (the x86 default NaN, 0xfff8000000000000); the C translation raises the flag the same way and stores
  that bit pattern explicitly, so the callout also returns it on AArch64.

## The C translation and the Highway version

As for asinh: `svml/svml_atanh_d_ha.c` (`svml_atanh8_ha`, `svml_atanh_ha_array`, `svml_datanh_ha_cout_rare`)
is bit-identical to `__svml_atanh8_ha` (`test_atanh8`) and has the same arithmetic instruction mix
(`make atanh-insns`). `hwy/hwy_atanh.cc` performs the same operations in the same order with the exact
emulations of `vrcp14pd` (twice), `vgetexppd` and the four 16-entry lookups; it is **bit-identical to
Intel's code** on every target with FMA (AVX2 and AVX3* under Intel SDE against the translation; NEON, SVE at
128–2048 bits and SVE2 under QEMU against x86-64 with `make atanh-xcheck`). Targets without FMA use
`std::atanh`. Deviations: none in the results; the Highway version may raise floating-point exception flags
that Intel's `{sae}` code suppresses.

Replacing the exact `vrcp14pd` by another 14-bit approximation (1/x truncated to 17 bits) changes 1 of
2,007,128 results: the rounded Rcp differs where rcp14 is within 2⁻¹⁴ of a rounding midpoint, and the other
(equally valid) reduction rounds differently.

## Tests (`make atanh-test`, `make atanh-avx512-test`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_datanh_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,007,122 inputs, 0 mismatches |
| `test_atanh8` | Intel's `__svml_atanh8_ha`: all 8 lanes bit-exact | 20,007,128 inputs, 0 mismatches (Intel SDE, `-spr`); needs AVX-512 (else skipped) |
| `hwy_test` | every Highway target: max error < 1 ulp vs `atanhl`, special results (±inf at ±1, NaN beyond, ±0), agreement between targets, bit-identity with the C translation when AVX-512 is available | all checks pass (SDE: AVX3, AVX3_DL, AVX3_SPR, AVX2 bit-identical to the translation) |
| `make atanh-insns` | arithmetic instruction mix of the compiled vector path vs the assembly | identical |

The emulated instructions themselves are tested by asinh's `test_approx14`.

Inputs (`tests/inputs.h`): ±0, ±1, ±inf, NaNs, denormals, ±DBL_MAX, 1 ± 1 ulp; ±100 ulps around ±1, ±0.5,
±2⁻⁵³, ±2⁻⁵⁴, ±0.25, ±0.75, ±2⁻⁵; |x| = j/32 and 1 − 2⁻ᵉ·(1 + j/32) with neighbours (boundaries of the
5-bit reciprocals and the tables of 1 + y and 1 − y); random values uniform in [−1, 1], with 1 − |x|
log-uniform in [2⁻⁵³, 1], log-uniform magnitudes below 1, and random bit patterns. All are built with IEEE
operations and bit manipulation only (no libm), so they are identical on every architecture.
