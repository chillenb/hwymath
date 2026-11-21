# expm1 (double precision)

Source studied: [`svml_z0_expm1_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_expm1_d_ha.s)
(AVX-512, high accuracy; vector entry `__svml_expm18_ha`, scalar callout `__svml_dexpm1_ha_cout_rare_internal`).
Accuracy and speed: [Performance.md](Performance.md).

| directory | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_expm1_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py`) |
| `hwy/` | portable Google Highway version of the same algorithm, bit-identical on every target with FMA |
| `tests/` | tests, comparison implementations (Intel SVML, SLEEF), test inputs |
| `tools/` | `gen_svml_constants.py`, `xcheck.sh` (AArch64 under QEMU) |
| `expm1_d.h` | declarations of the C kernels |
| `function.mk` | build rules and the `expm1-*` make targets |

## Algorithm

Vector path (8 lanes). Every operation rounds to nearest; the only AVX-512-specific instructions are
`vscalefpd` and `vpermt2pd`.

- t = x·log2(e) + shifter (an FMA), with shifter = 1.5·2^48 + 1023. Then N = t − shifter is x·log2(e)
  rounded to a multiple of 1/16, and the low 4 bits of t give j.
- Th[j], Tlr[j] come from 16-entry tables (`vpermt2pd`), with 2^(j/16) = Th·(1 + Tlr).
- Th' = `vscalefpd(Th, max(N, −128))` = Th·2^floor(N), clamped so that it stays normal for very
  negative x.
- R = x − N·ln2_hi (one FMA; ln2_hi is the full double ln 2), then R' = R − N·ln2_lo, and Rl = (R − R')
  − N·ln2_lo, its rounding error. R' and Rl have their top exponent bit cleared (`and` with
  0xbfff…), so |R'| < 2 even for −inf and huge |x|.
- P = R'²·(c2 + c3 R' + … + c8 R'^6) + Tlr·R' + (Rl + Tlr), from Taylor-like coefficients.
- Th' − 1 is a Fast2Sum, where the operand order follows Th' < 1 (`vblendmpd`): tm1 = Th' − 1 and its
  error Tl1. Th'·R' is a double-double: hi = Th'·R' + tm1, lo = Th'·R' − (hi − tm1).
- expm1(x) = hi + ((P·Th' + Tl1) + lo), and the sign bit of x is ORed in, so expm1(−0) = −0.
- Lanes with x ≥ 707.703… (NLT_UQ: also +inf and NaN) take the scalar callout. Large negative x
  stays in the vector path and gives −1, including −inf.

The scalar callout is SVML's scalar exp (the same code and the same 64-entry `_imldExpHATab` table as
exp's callout), preceded by `x < 0 → −1`. It is only meant for the lanes that the vector path sends
it, where exp(x) − 1 rounds to exp(x). Its comparison `x < 0` is a signalling `comisd`, so a quiet
NaN raises the invalid flag. The C version keeps this, because it uses a C relational operator.

## The C translation (`svml/svml_expm1_d_ha.c`)

`svml_expm18_ha(__m512d)` translates the assembly into C with AVX-512F/DQ intrinsics. It keeps the
algorithm, the constants (bit-exact hex-float literals), the operation order, the operand order of
the non-commutative operations (`vmaxpd`, `vblendmpd`) and each instruction's rounding/SAE mode, so its
results are bit-identical to Intel's (`test_expm18`). `vandpd`/`vorpd` need AVX-512DQ. An empty
`asm` barrier on the sign mask keeps GCC from fusing the `and` and the final `or` into one
`vpternlogq`. That fusion gives the same results, but the instruction mix would no longer match the
assembly's in `make expm1-insns`. `svml_dexpm1_ha_cout_rare()` is the scalar callout in plain C, and
it is bit-identical to the assembly's, including return codes and FP exception flags (`test_rare`).
Its scalar table is a local symbol named like exp's table, so each function's assembly stays in its
own object.

## Highway version (`hwy/hwy_expm1.cc`)

The same algorithm with [Google Highway](https://github.com/google/highway) 1.4 and run-time dispatch.
On targets with true FMA (`HWY_NATIVE_FMA`: AVX2, AVX3*, NEON, SVE/SVE2, …) the results are
**bit-identical to Intel's AVX-512 code**. The AVX-512-only operations are emulated exactly:

| AVX-512 feature | emulation | why it is exact |
|---|---|---|
| `vscalefpd(Th, max(N, −128))` | Th · 2^floor(max(N, −128)), where the power of two is built from the bits of max(N, −128) + shifter (`(bits >> 4) << 52`: the shifter's mantissa holds 16·(N + 1023)) | max(N, −128) is a multiple of 1/16 in [−128, 1022) for the vector-path lanes, so the addition is exact. Th·2^M is normal, so the product is exact, like `vscalefpd` |
| `vpermt2pd` (16-entry tables) | ≥ 16 lanes: `TableLookupLanes` on the table repeated to 64 entries; 8 lanes: `TwoTablesLookupLanes`; 4 lanes: two `TwoTablesLookupLanes` + select on bit 3; 2 lanes: four + two selects. No gathers | table lookups |
| `vmaxpd {sae}`, `vblendmpd`, `vcmppd`, `vandpd`/`vorpd` | `Max`, `IfThenElse`, `Lt`, `And`/`Or` on the bits | N is never NaN in the vector-path lanes, so the NaN rules of `vmaxpd` do not matter |

Everything else is the same sequence of round-to-nearest FMAs (`MulAdd`, `NegMulAdd`, `MulSub`),
multiplications and additions. Lanes with !(x < 707.703…) use the C callout from `svml/` (an
out-of-line `cold` function). Targets without FMA (SSE2/SSSE3/SSE4) use `std::expm1`.

**AVX-512 targets** (AVX3, AVX3_DL, AVX3_ZEN4, AVX3_SPR, AVX10_2; `HWYMATH_NATIVE_AVX512` from
`common/include/hwymath/hwy/avx512_native-inl.h`) use the instructions instead of the emulations above:
`vpermt2pd` directly on the bits of t (the instruction ignores all but the low 4 index bits, so the
`And(t, 15)` is gone), `vscalefpd(Th, vmaxpd(N, −128))`, and the callout lanes are found with one
`vcmppd` NLT_UQ (`_mm512_cmp_pd_mask`; `Not(Lt(...))` compiled to a compare and a `knotb`). The raw
`vpermt2pd` and the NLT_UQ compare are small local wrappers in `hwy_expm1.cc` (Highway's lookups require
in-range indices). The AVX3 loop now has 37 vector arithmetic/mask instructions per 8 values (before: 42),
Intel's main path 38 (its `vandpd` + `vorpd` of the sign are one `vpternlogq` here, its `kmovw` + `testl`
one `kortestb`); the remaining difference is register copies (Highway keeps the constants in registers,
Intel reloads them).

Default NaN: for x = −inf (and for x = −DBL_MAX, where x·log2(e) overflows), R or Rl is an invalid
operation (inf − inf). x86 produces the default NaN 0xfff8…, and AArch64 produces 0x7ff8…. After
the `0xbfff…` mask, R' or Rl is then −1.5 on x86 and +1.5 on AArch64. These values are scaled by
Th' = 2^−128 or less, so both architectures still return exactly −1. `make expm1-xcheck` confirms
bit-identity, and these inputs are in `tests/inputs.h`.

## Tests (`make expm1-test`, `make expm1-avx512-test`, `make expm1-insns`, `make expm1-xcheck`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dexpm1_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,014,924 inputs, 0 mismatches |
| `test_expm18` | Intel's `__svml_expm18_ha`: all 8 lanes bit-exact | 2,014,928 inputs, 0 mismatches (Intel SDE `-spr`) |
| `hwy_test` | every Highway target: error < 1 ulp, special values, identical results on all FMA targets, and with AVX-512 bit-identity with the C translation | AVX3_SPR, AVX3_DL, AVX3, AVX2 bit-identical to the C translation (SDE, 2,014,924 inputs) |
| `make expm1-insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector path vs the assembly | identical |
| `make expm1-xcheck` | static AArch64 NEON, SVE and SVE2 builds under QEMU (SVE at 128/256/512/2048 bits: 2/4/8/32 lanes, so every lookup variant) vs x86-64 AVX2 | bit-identical (2,014,924 inputs) |

Inputs (`tests/inputs.h`): special values (signed zeros, infinities, quiet and signalling NaNs with
payloads, extremes such as ±DBL_MAX and −1e100), ±50 ulps around every threshold (707.703, 709.78,
the N = −128 clamp at about −88.72, the −1 saturation at about −37.43, ±ln2/32 where the Fast2Sum
order switches, …), ±2^e and neighbours for every e from −1074 to 10, table and rounding boundaries
(16M + j/2)·ln2/16 and their neighbours, and random values (uniform in [−750, 750] and [−1, 1], and
random bit patterns).

`hwy_test --hash FILE` writes the result hashes of every FMA target and, when AVX-512 is available,
of the C translation (`C-AVX512`). All of them share the variant `svml`, so the CI cross-check
compares the translation and every Highway target on every machine.

## Comparisons

`make expm1-hwy` runs the Highway targets against Intel's libsvml (`__svml_expm14_ha_l9`,
`__svml_expm18_ha_z0`; conda-forge `intel-cmplr-lib-rt`, closed source, used for comparison only),
SLEEF (u10), glibc's scalar `expm1` and Highway's own `Expm1`. glibc 2.28 (the conda sysroot) has no
vector expm1.
