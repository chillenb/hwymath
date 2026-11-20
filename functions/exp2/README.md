# exp2 (double precision)

Source studied: [`svml_z0_exp2_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_exp2_d_ha.s)
(AVX-512, high accuracy; vector entry `__svml_exp28_ha`, scalar callout `__svml_dexp2_ha_cout_rare_internal`).
Accuracy and speed: [Performance.md](Performance.md).

| directory | contents |
|---|---|
| `svml/` | bit-identical C/AVX-512 translation of the assembly, with its constants (`svml_exp2_d_ha_constants.inc`, generated from the assembly by `tools/gen_svml_constants.py`) |
| `hwy/` | portable Google Highway version of the same algorithm, bit-identical on every target with FMA |
| `tests/` | tests, comparison implementations (Intel SVML, SLEEF), test inputs |
| `tools/` | `gen_svml_constants.py`, `xcheck.sh` (AArch64 under QEMU) |
| `exp2_d.h` | declarations of the C kernels |
| `function.mk` | build rules and the `exp2-*` make targets |

## Algorithm

Vector path (8 lanes), with j = floor(16x) mod 16:

- `r = vreducepd(x, 0x41) = x − floor(16x)/16` ∈ [0, 1/16). The immediate selects 4 fraction bits and
  rounding toward −∞, for the rounding to an integer and for the subtraction. The subtraction is exact
  except for −1/16 < x < 0, where r = x + 1/16 is rounded down.
- `t = x + 1.5·2^48` with `{rd-sae}` (round toward −∞). The low 4 mantissa bits of t are j (`vpandq` with 15).
- 2^(j/16) = Th[j] + Tl[j], from two 16-entry tables (`vpermt2pd`, two sources of 8).
- y = Th + (Tl + (Th·r)·P(r)), where P(r) = c1 + c2 r + … + c6 r^5 (Horner, FMAs).
- `exp2(x) = vscalefpd(y, x) = y · 2^floor(x)`.
- Lanes with |x| ≥ 1021.99951171875 (and ±inf) take the scalar callout. NaN compares unordered and
  never takes the callout: the vector path returns the quieted NaN.

Scalar callout (64-entry {hi, lo} table, x = M + j/64 + r with |r| ≤ 1/128, degree-5 polynomial):
- inf/NaN: exp2(−inf) = +0, otherwise x·x;
- x ≥ 1024: overflow (DBL_MAX², code 3);
- x ≤ −1075: underflow (tiny², code 4);
- x < −1022 (denormal results): it scales by 2^(M+100) and then by 2^−100, which is a double
  rounding (code 4);
- otherwise: 2^M · p, with M = 1024 split as 2^1023 · 2.

Unlike exp's callout, this one has no extra-precision path for denormal results. That is where the
largest error occurs (0.752 ulp at x ≈ −1022.44).

## The C translation (`svml/svml_exp2_d_ha.c`)

`svml_exp28_ha(__m512d)` translates the assembly into C with AVX-512F/DQ intrinsics. It keeps the
algorithm, the constants (bit-exact hex-float literals), the operation order and each instruction's
rounding mode (`_mm512_reduce_round_pd(x, 0x41, …)`, `_mm512_add_round_pd(…, _MM_FROUND_TO_NEG_INF)`),
so its results are bit-identical to Intel's (`test_exp28`). `svml_exp2_ha_array(x, y, n)` is the
array loop, with a masked tail. `svml_dexp2_ha_cout_rare()` is the scalar callout in plain C (no FMA
contraction, `volatile` where flags must be raised at run time). It is bit-identical to the
assembly's callout, including return codes and FP exception flags (`test_rare`). Only the vector
functions are compiled for AVX-512 (`target("avx512f,avx512dq")`). `vreducepd` and `vandpd` need
AVX-512DQ.

## Highway version (`hwy/hwy_exp2.cc`)

The same algorithm with [Google Highway](https://github.com/google/highway) 1.4 and run-time dispatch.
On targets with true FMA (`HWY_NATIVE_FMA`: AVX2, AVX3*, NEON, SVE/SVE2, …) the results are
**bit-identical to Intel's AVX-512 code**. Each AVX-512-only operation is emulated exactly:

| AVX-512 feature | emulation | why it is exact |
|---|---|---|
| `vreducepd {sae}`, imm 0x41 (round down) | f = Floor(16x)/16; r = x − f (RN); err = x − (r + f); if err < 0, step r one ulp down (integer −1 on its bits) | f ≤ x. For x ≥ 0 or x ≤ −1/16 the subtraction is exact, and so err = 0. For −1/16 < x < 0, \|f\| > \|x\|, so Fast2Sum gives the exact error, and RD = RN − 1 ulp exactly when RN rounded up (r ≈ 1/16 > 0) |
| `vaddpd {rd-sae}` x + 1.5·2^48 (only its low 4 bits are used) | j = (Floor(16x) + 1.5·2^52 + 16·1023) & 15 | RD(x + 1.5·2^48) = 1.5·2^48 + Floor(16x)/16 for \|x\| < 2^47 |
| `vscalefpd(y, x)` = y·2^floor(x) | multiply by 2^floor(x), whose exponent field is bits 4..15 of the same shifted Floor(16x) | for the vector-path lanes (\|x\| < 1022), 2^floor(x) and y·2^floor(x) are normal, so the product is exact, like `vscalefpd` |
| `vpermt2pd` (16-entry tables) | ≥ 16 lanes: `TableLookupLanes` on the table repeated to 64 entries; 8 lanes: `TwoTablesLookupLanes`; 4 lanes: two `TwoTablesLookupLanes` + select on bit 3; 2 lanes: four + two selects. No gathers | table lookups |
| `vcmppd` GE_OQ on \|x\| | `!(Abs(x) < thr)` → scalar fix-up, which also catches NaN | for NaN, the callout returns x·x, the same quieted NaN as Intel's vector path |

Lanes outside the vector range use the C callout from `svml/` (an out-of-line `cold` function). Targets
without FMA (SSE2/SSSE3/SSE4) use `std::exp2`.

Notes:
- −0 versus +0: for multiples of 1/16, r is zero, and its sign may differ between `vreducepd`
  (x − x rounded down is −0) and the emulation (+0). The sign does not reach the result: Th·r·P(r) +
  Tl is a sum with a zero, and Tl ≠ 0 except for j = 0, where Tl = +0 and (±0) + (+0) = +0. The
  bit-identity tests include these inputs.
- The tables are `constexpr` (the generated `.inc` uses an `SVML_TABLE` storage-class macro that
  `hwy_exp2_constants.h` sets to `static constexpr`). Run-time-initialized padded tables would put
  init guards into the hot loop.
- AVX2 cost: the two 16-entry lookups take 8 `vpermd` and 6 blends per vector. This is the main
  reason why the Highway AVX2 version is slower than Intel's own AVX2 exp2, which is a different code
  path (Performance.md).

## Tests (`make exp2-test`, `make exp2-avx512-test`, `make exp2-insns`, `make exp2-xcheck`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dexp2_ha_cout_rare_internal`: results, return codes, FP exception flags | 20,016,321 inputs, 0 mismatches |
| `test_exp28` | Intel's `__svml_exp28_ha`: all 8 lanes bit-exact | 2,016,328 inputs, 0 mismatches (Intel SDE `-spr`) |
| `hwy_test` | every Highway target: error < 1 ulp, special values, identical results on all FMA targets, and with AVX-512 bit-identity with the C translation | AVX3_SPR, AVX3_DL, AVX3, AVX2 bit-identical to the C translation (SDE, 2,016,321 inputs) |
| `make exp2-insns` | arithmetic instruction mix (mnemonic + rounding mode) of the compiled vector path vs the assembly | identical |
| `make exp2-xcheck` | static AArch64 NEON, SVE and SVE2 builds under QEMU (SVE at 128/256/512/2048 bits: 2/4/8/32 lanes, so every lookup variant) vs x86-64 AVX2 | bit-identical (2,016,321 inputs) |

Inputs (`tests/inputs.h`): special values (signed zeros, infinities, quiet and signalling NaNs with
payloads, extremes), ±50 ulps around every threshold of both algorithms (±1021.9995, 1024, 1023.992,
−1022, −1074, −1075, ±1/16, …), ±2^e for every e from −1074 to 0 (the rounded `vreducepd` subtraction),
table boundaries M + j/64 and their neighbours, and random values (uniform in [−1100, 1100] and
[−1, 1], and random bit patterns). The assembly is read from the submodule. A copy with the scalar
routine made global is assembled into `build/exp2/`. Disabling the round-down correction of the
`vreducepd` emulation makes `hwy_test` fail (25 of 200k inputs differ), so the inputs cover it.

`hwy_test --hash FILE` writes the result hashes of every FMA target and, when AVX-512 is available,
of the C translation (`C-AVX512`). All of them share the variant `svml`, so the CI cross-check
compares the translation and every Highway target on every machine.

## Comparisons

`make exp2-hwy` runs the Highway targets against Intel's libsvml (`__svml_exp24_ha_l9`, `__svml_exp28_ha_z0`;
conda-forge `intel-cmplr-lib-rt`, closed source, used for comparison only), SLEEF (u10), glibc's scalar
`exp2` and Highway's own `Exp2`. glibc 2.28 (the conda sysroot) has no vector exp2. libsvml's
AVX-512 exp2 is not the same code as the published assembly: its results differ (max 0.763 ulp versus
0.748 ulp on the same 2M inputs).
