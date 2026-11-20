# cos (double precision)

Source studied: [`svml_z0_cos_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_cos_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md).

cos is sin's twin: same data layout and tables, same large-argument path. This directory holds what is
specific to cos; the shared parts (large-argument reduction in C and Highway, test templates and inputs,
constant generator, `trig_common.mk`) live in [`functions/sin`](../sin/README.md), whose README also
describes the algorithm, the emulation of each AVX-512 feature and the test inputs in detail.

| path | contents |
|---|---|
| `svml/svml_cos_d_ha.c` | bit-identical C/AVX-512 translation of the assembly and its scalar callout |
| `svml/svml_cos_d_ha_constants.inc` | its constants (`functions/sin/tools/gen_svml_constants.py cos`) |
| `hwy/hwy_cos.cc` | Highway version |
| `tests/` | `test_rare.c`, `test_cos8.c`, `hwy_test.cc`, `extras.c` (thin wrappers of the templates in `functions/sin/tests`) |
| `cos_d.h` | declarations of the C translation |

## Algorithm

**Main path** (|x| <= 2^24): N = round(x*16/pi) (`y = x*16/pi + 1.5*2^52`), R + Rlo = x - N*pi/16 in
three parts, and with B = N*pi/16 and the tables of sin (Sh + Sl = sin B, Ch + Cl = cos B):
cos(B + R) = Ch - R*Sh [with its exact rounding error] + R^2*(Ch*Cp - R*Sh*Sp) + Cl
&minus; Rlo*(Sh + Ch*R) &minus; R*Sl, with cos R - 1 = R^2*Cp(R^2), sin R - R = R^3*Sp(R^2); bit 4 of N flips the sign. No special
case for -0 (cos is even).

**Large path** (|x| > 2^24, +-inf, NaN): the Payne-Hanek reduction of sin, applied to x itself (with its
sign) instead of |x|, and the same evaluation with cos's 512-entry table: its rows {T0, T1, T2, T3} hold
T2 + T3 = cos B and T0 + T1 = -sin B (T1 a power of 2), so that the same formula
f(B + r) = T2 + (T0 + T1)*r + T2*(cos r - 1) + (T0 + T1)*(sin r - r) + ... gives cos(B + r); no final
sign change.

**Callout**: cos(+-inf) = NaN (invalid, return code 1), cos(NaN) = x*x; same routine as sin's.

## The C translation and the Highway version

As for sin (see [its README](../sin/README.md)): `svml_cos8_ha(__m512d)`, `svml_cos_ha_array()`,
`svml_dcos_ha_cout_rare()`; the Highway version emulates every AVX-512 feature exactly (16-entry
`vpermt2pd` lookups by permutes, `vpmullq` by `MulEven`, the gathers of the large path by scalar loads),
so it is bit-identical to Intel's AVX-512 code on every FMA target; targets without FMA use `std::cos`.
**Deviations: none** on FMA targets; `make cos-xcheck` confirms bit-identical results on AArch64
NEON/SVE/SVE2 under QEMU (SVE at 128 to 2048 bits). Intel's libsvml (`__svml_cos8_ha_z0`) agrees with the translation on
all test inputs except cos(+-inf) (+NaN instead of -NaN).

## Tests (`make cos-test`, `make cos-avx512-test RUN="$SDE -spr --"`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's `__svml_dcos_ha_cout_rare_internal`: results, whether written, return codes, FP flags | 20,035,562 inputs, 0 mismatches |
| `test_cos8` | Intel's `__svml_cos8_ha`, bit for bit (Intel SDE) | 20,035,568 inputs (7,350,780 in the large path), 0 mismatches |
| `hwy_test` | every Highway target: error < 1 ulp, special values, agreement between targets, bit-identity with the C translation (SDE: AVX3_SPR, AVX3_DL, AVX3, AVX2) | all checks pass |
| `make cos-insns` | instruction mix of the main path vs the assembly | identical |

Inputs: `functions/sin/tests/inputs.h` (described in sin's README).
