# exp (double precision): performance

Accuracy and throughput of the exp implementations in this directory and of the libraries they are
compared with, on three machines. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `expl` (`long double`), over the
  test inputs of `tests/inputs.h` (special values, every threshold of the algorithms with neighbouring
  doubles, table boundaries, random values in [-750, 750] and random bit patterns); statistics over the
  inputs with a finite, nonzero result. "not CR": results that are not correctly rounded (error > 0.5 ulp).
  The x86-64 `long double` has a 64-bit mantissa, which misjudges some results lying almost exactly
  halfway between two doubles, so the x86-64 "not CR" percentages are slightly optimistic; AArch64's
  `long double` is binary128 (113-bit mantissa).
- **Throughput**: nanoseconds per value for 4096 random doubles in [-700, 700] (cache-resident), best of
  7 runs of 2000 repetitions. Array functions process the whole array in one call (constants hoisted out
  of the loop); per-call numbers call a 4- or 8-lane function once per vector, the way vectorizing
  compilers and NumPy call SVML. Intel's libsvml can only be measured per call.
- **Toolchain**: conda-forge GCC 15.3 via pixi. It links against a glibc 2.28 sysroot, so the scalar
  `exp` baseline binds to the older, slower `exp@GLIBC_2.2.5` symbol version; scalar numbers are
  pessimistic. CI machines are shared virtual machines: expect some run-to-run variation in the timings
  (the accuracy figures are exact).

Programs: `make exp-hwy` (`tests/hwy_test.cc`: Highway targets and every other implementation available),
`make exp-compare` / `make exp-bench` (x86-64 AVX2 variants, the latter also per call).

## Summary

| | Intel Xeon 6975P-C (AVX-512) | | Neoverse-N2 (SVE2, 128-bit) | | Core i5-9500 (AVX2) | |
|---|---|---|---|---|---|---|
| | max ulp | ns/value | max ulp | ns/value | max ulp | ns/value |
| **Highway hp, best target** | **0.505** (AVX3, 8-entry) | **0.55** | **0.511** (SVE2) | **3.69** | **0.511** (AVX2) | **1.07** |
| Intel SVML HA, widest | 0.800 (AVX-512) | 0.40 | — | — | 0.521 (AVX2) | 1.03 |
| SLEEF u10, widest | 0.877 (AVX-512) | 0.69 | 0.877 (SVE) | 4.33 | 0.877 (AVX2) | 1.50 |
| glibc libmvec, widest | 3.022 (AVX-512) | 0.43 | — | — | 3.008 (AVX2) | 2.32 |
| Highway contrib `Exp`, best target | 1.030 | 0.65 | 1.030 | 4.11 | 1.013 | 1.35 |

## Intel Xeon 6975P-C (Granite Rapids), GitHub Actions `ubuntu-24.04`

CI run [44747421851](https://github.com/chillenb/hwymath/actions/runs/44747421851) (commit `113a127`),
20,002,324 inputs, `long double` with 64-bit mantissa. AVX-512 F/CD/BW/DQ/VL/IFMA/VBMI/VBMI2/VNNI/BITALG/
VPOPCNTDQ/BF16/FP16.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|---|
| Highway AVX3_SPR | hp, 8-entry | 8 | 0.5054 | 0.0081 % | 0 | 0.57 |
| Highway AVX3_DL | hp, 8-entry | 8 | 0.5054 | 0.0081 % | 0 | 0.57 |
| Highway AVX3 | hp, 8-entry | 8 | 0.5054 | 0.0081 % | 0 | 0.55 |
| Highway AVX2 | hp, 4-entry | 4 | 0.5107 | 0.0557 % | 0 | 0.91 |
| Highway SSE4 / SSSE3 / SSE2 | `std::exp` | 2 | 0.5063 | 0.0362 % | 0 | 5.3–5.5 |
| Highway contrib `Exp`, AVX3 | | 8 | 1.0303 | 6.1069 % | 3 | 0.65 |
| Highway contrib `Exp`, AVX2 | | 4 | 1.0303 | 6.1069 % | 3 | 1.13 |
| Highway contrib `Exp`, SSE4 | | 2 | 0.9805 | 6.3668 % | 0 | 3.08 |
| C hp AVX2 (array) | | 4 | 0.5107 | 0.0557 % | 0 | 0.87 |
| C fast AVX2 (array) | | 4 | 0.6455 | 1.0915 % | 0 | 0.85 |
| C bit-identical SVML port AVX2 (array) | | 4 | 0.7217 | 3.0555 % | 0 | 1.71 |
| C SVML translation AVX-512 (array) | | 8 | 0.7217 | 3.0555 % | 0 | 0.49 |
| Intel libsvml HA AVX2 (`__svml_exp4_ha_l9`) | | 4 | 0.5225 | 0.3213 % | 0 | 0.84 |
| Intel libsvml HA AVX-512 (`__svml_exp8_ha_z0`) | | 8 | 0.8000 | 3.0768 % | 0 | 0.40 |
| SLEEF u10 AVX2 | | 4 | 0.8770 | 3.9164 % | 0 | 1.30 |
| SLEEF u10 AVX-512 | | 8 | 0.8770 | 3.9164 % | 0 | 0.69 |
| glibc libmvec AVX2 | | 4 | 3.0220 | 28.13 % | 1,104,552 | 0.67 |
| glibc libmvec AVX-512 | | 8 | 3.0220 | 28.13 % | 1,104,552 | 0.43 |
| glibc `exp` (scalar) | | 1 | 0.5063 | 0.0362 % | 0 | 5.31 |

Per call (`make exp-bench`): Intel AVX2 HA 0.84, Intel AVX2 LA 0.73, C hp 0.94, C fast 0.90, bit-identical
port 1.77 ns/value.

`test_exp8` on this machine: the AVX-512 C translation and the AVX2 port are both bit-identical to
Intel's `__svml_exp8_ha` assembly on 100,002,328 inputs.

## Neoverse-N2, GitHub Actions `ubuntu-24.04-arm`

Same CI run, 20,002,324 inputs, `long double` binary128. SVE/SVE2 with 128-bit vectors (2 doubles).

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|---|
| Highway SVE2_128 | hp, 4-entry | 2 | 0.5109 | 0.0694 % | 0 | 3.79 |
| Highway SVE2 | hp, 4-entry | 2 | 0.5109 | 0.0694 % | 0 | 3.69 |
| Highway SVE | hp, 4-entry | 2 | 0.5109 | 0.0694 % | 0 | 4.34 |
| Highway NEON | hp, 4-entry | 2 | 0.5109 | 0.0694 % | 0 | 4.27 |
| Highway NEON_WITHOUT_AES | hp, 4-entry | 2 | 0.5109 | 0.0694 % | 0 | 4.27 |
| Highway contrib `Exp` (SVE2) | | 2 | 1.0300 | 6.1256 % | 3 | 4.14 |
| Highway contrib `Exp` (NEON) | | 2 | 1.0300 | 6.1256 % | 3 | 4.61 |
| SLEEF u10 SVE | | 2 | 0.8773 | 3.9349 % | 0 | 4.33 |
| SLEEF u10 AdvSIMD | | 2 | 0.8773 | 3.9349 % | 0 | 4.59 |
| glibc `exp` (scalar) | | 1 | 0.5061 | 0.0508 % | 0 | 7.01 |

## Intel Core i5-9500 (Coffee Lake, AVX2, no AVX-512)

Local runs; the "Gather Data Sampling" microcode mitigation is active on this CPU, which makes AVX2
gathers very slow (none of the final kernels use gathers).

AVX2 variants (`make exp-compare`, `make exp-bench`; 20,002,324 inputs, 64-bit-mantissa `long double`):

| implementation | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|
| bit-identical SVML port (AVX2) | 0.722 | 3.05 % | 0 | 1.94 (1.69 array) |
| fast variant, K=3 (default) | 0.635 | 1.09 % | 0 | 1.10 (0.89 array) |
| fast variant, K=2, degree 8 | 0.789 | 2.04 % | 0 | (0.78 array, before the cold-callout change) |
| **hp variant, K=2, Q to r^9 (default)** | **0.511** | **0.055 %** | 0 | 1.25 (1.04 array) |
| hp variant, K=3, Q to r^8 | 0.5068 | 0.0079 % | 0 | 1.54 (1.17 array) |
| Intel `__svml_exp4_ha_l9` (AVX2 HA) | 0.521 | 0.32 % | 0 | 1.03 |
| Intel `__svml_exp4_l9` (AVX2 LA) | 3.008 | 28.1 % | 1,104,775 | 0.84 |
| glibc `exp` (scalar) | 0.507 | 0.036 % | 0 | 4.98 |
| glibc libmvec `_ZGVdN4v_exp` (AVX2) | 3.008 | 28.1 % | 1,104,775 | 2.32 |
| SLEEF 3.9.0 `Sleef_expd4_u10avx2` | 0.877 | 3.91 % | 0 | 1.50 |

Highway (`make exp-hwy`, `make exp-xcheck`):

| target | algorithm | max ulp | not CR | ns/value (array) | check |
|---|---|---|---|---|---|
| AVX2 | hp, 4-entry | 0.511 | 0.055 % | 1.07 (C: 1.03) | bit-identical to the C hp variant (20,002,324 inputs) |
| SSE4 / SSSE3 | `std::exp` | 0.507 | 0.036 % | 4.75 | |
| AArch64 NEON (QEMU) | hp, 4-entry | 0.510* | 0.068 %* | | bit-identical to x86-64 AVX2 |
| AArch64 SVE, SVE2 at 128/256/512/2048 bits (QEMU) | hp, 4-entry | 0.510* | 0.068 %* | | bit-identical to x86-64 AVX2 |
| Highway contrib `Exp` (AVX2) | | 1.013 | 6.09 % | 1.35 | |

Highway port of Intel's algorithm (`hwy/hwy_exp_svml.cc`; `make exp-hwy`, 20,002,324 inputs, and
`make exp-xcheck`). Timings from the i5-9500 only (AVX2; there is no AVX-512 on this machine, so the
AVX-512 targets were only checked for bit-identity, under Intel SDE `-spr`, where timings are
meaningless; CI numbers are to be added at integration):

| target | algorithm | max ulp | not CR | ns/value (array) | check |
|---|---|---|---|---|---|
| AVX2 | svml | 0.7217 | 3.0555 % | 1.56 | bit-identical to the C AVX2 port (1.70 ns/value in the same run) |
| AVX3_SPR / AVX3_DL / AVX3 (SDE) | svml | 0.6987* | 3.0736 %* | | bit-identical to the C/AVX-512 translation and to Intel's `__svml_exp8_ha` |
| SSE4 / SSSE3 / SSE2 | `std::exp` | 0.5063 | 0.0362 % | 6.8–7.0 | |
| AArch64 NEON, SVE, SVE2 at 128/256/512/2048 bits (QEMU) | svml | 0.6989† | 3.094 %† | | bit-identical to x86-64 AVX2 |

\* 2,002,324 inputs. † 2,002,324 inputs, binary128 reference.

Same run, for comparison: Intel `__svml_exp4_ha_l9` 1.04 ns/value (0.5225 ulp), SLEEF u10 AVX2 1.52
(0.877 ulp), Highway hp AVX2 1.15 (0.5107 ulp). The Highway port is ~8 % faster than the C AVX2 port
(fewer instructions: no `and` of r, no masking of the table index on AVX2, 2^M from two shifts).

\* 2,002,324 inputs, binary128 reference; the same results give 0.5098 ulp / 0.054 % with the x86-64
reference.

## Observations

- The Highway hp exp is the most accurate vectorized exp measured here on every machine (max ≤ 0.511 ulp,
  not correctly rounded for < 0.07 % of inputs; 0.505 ulp / 0.008 % with the 8-entry AVX-512 table), and its
  results are bit-identical across x86-64 AVX2, AArch64 NEON and SVE/SVE2 (CI cross-check job).
- Speed: on par with Intel's AVX2 HA exp in array form; on AVX-512 about 35 % slower than Intel's AVX-512
  HA exp, which is much less accurate (0.800 ulp). On Neoverse-N2 it is the fastest of the vectorized exps
  measured, ahead of SLEEF's SVE and AdvSIMD versions.
- Intel's current libsvml uses different exp algorithms than the open-source numpy/SVML AVX-512 code: its
  AVX2 HA version is more accurate (0.52 ulp) and its AVX-512 HA version less accurate (0.80 ulp) than
  the open-source AVX-512 HA algorithm (0.72 ulp), which this repository's bit-identical port reproduces.
- glibc's libmvec AVX2 exp has exactly the error profile of Intel's AVX2 LA exp (it was contributed by Intel).
- Highway's `hwy/contrib/math` `Exp` exceeds its documented 1 ulp on 3 of the inputs (max 1.030 ulp,
  x86-64 and AArch64 alike).
- The remaining speed gap of the bit-identical AVX2 port is the cost of staying identical to the AVX-512
  algorithm (16-entry table lookups, round-toward-zero emulation); its first version used gathers and ran
  at 4.26 ns/value on the i5-9500 because of the GDS microcode mitigation.
