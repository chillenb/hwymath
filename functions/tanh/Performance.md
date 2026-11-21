# tanh (double precision): performance

Accuracy and throughput of the tanh implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result against `tanhl` (`long double`), over the
  inputs of `tests/inputs.h` (special values, ±50 doubles around every subinterval boundary, the
  saturation points and the callout threshold, uniform values in [−30, 30] and [−2, 2], log-uniform
  magnitudes 2^−60 ... 2^6, random bit patterns); statistics over the inputs with a finite, nonzero
  reference, and exact agreement required for NaN/inf/zero results ("special" column: mismatches). "not
  CR": error > 0.5 ulp. The x86-64 `long double` has a 64-bit mantissa, which misjudges some results lying
  almost exactly halfway between two doubles, so the x86-64 "not CR" percentages are slightly optimistic;
  the AArch64 run (binary128 reference) shows the difference.
- **Throughput** (`make tanh-hwy`): nanoseconds per value for 4096 random doubles uniform in [−20, 20]
  (cache-resident, no callout lanes), best of 7 runs of 2000 repetitions, on a Core i5-9500 (Coffee Lake,
  AVX2, no AVX-512) shared with other jobs: expect a few percent of noise. Intel's libsvml and SLEEF are
  called per vector of 4 lanes through the same array loop. SDE timings are meaningless and not given;
  AVX-512 timings will come from CI.
- **Toolchain**: conda-forge GCC 15.3, Highway 1.4.0, SLEEF 3.9.0, Intel libsvml from
  `intel-cmplr-lib-rt` 2026.2.0 (comparison only). glibc's scalar `tanh` binds to the glibc 2.28 sysroot's symbol.

## Core i5-9500 (Coffee Lake, AVX2), local

`make tanh-hwy`: 20,005,683 inputs, x86-64 `long double` reference.

| implementation | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm) | 4 | **0.6436** | **0.7872 %** | 0 | **3.17** |
| Highway SSE4 / SSSE3 / SSE2 (`std::tanh`) | 2 | 2.1421 | 12.6490 % | 490,299 | 16.1–16.6 |
| Intel libsvml HA AVX2 (`__svml_tanh4_ha_l9`) | 4 | 0.7539 | 0.6562 % | 0 | 3.35 |
| SLEEF u10 AVX2 (`Sleef_tanhd4_u10avx2`) | 4 | 1.0005 | 0.3542 % | 1,841 | 5.59 |
| glibc `tanh` (scalar) | 1 | 2.1421 | 12.6490 % | 490,299 | 16.04 |

All Highway targets with FMA give bit-identical results, so the accuracy figures of the Highway
AVX2 row hold for every such target (AVX-512, NEON, SVE).

## AVX-512 (Intel SDE, `-spr`), accuracy only

`make tanh-avx512-test`: 2,005,683 inputs.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C translation AVX-512 = Highway AVX3_SPR/AVX3_DL/AVX3/AVX2 (bit-identical) | 0.6357 | 0.7819 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_tanh8_ha_z0`) | 0.6357 | 0.7819 % | 0 |
| SLEEF u10 AVX-512 | 1.0005 | 0.3628 % | 181 |

## AArch64 (QEMU), binary128 reference

`common/tools/xcheck.sh` (as in `make tanh-xcheck`, with `NR` = 300,000): 300,000 random inputs plus the special
and threshold inputs: NEON, SVE and SVE2 at 128 to 2048 bits are bit-identical to x86-64; against the
binary128 `long double` reference the error is max 0.6340 ulp, 0.8163 % not correctly rounded (the same
results measured with the x86-64 reference: 0.6343 ulp, 0.8098 %).

## Observations

- Highway's AVX2 version is on par with (about 5 % faster than) Intel's AVX2 HA tanh, which uses a
  different algorithm there (max 0.754 ulp); Intel's AVX-512 libsvml tanh has exactly the statistics of the
  open-source algorithm ported here. The 19 table lookups decide the speed: 5.92–5.95 ns/value with the
  permute/blend tree, 3.25–3.36 with the row lookups; the remaining time is mostly the 16-FMA Horner chain.
- Max error 0.64 ulp; 0.79 % not correctly rounded. Saturation (±1 from |x| ≈ 19.06), tanh(±0) = ±0,
  tanh(±inf) = ±1 and tiny arguments are exact.
- glibc's scalar `tanh` reaches 2.14 ulp (490,299 results above 1 ulp); SLEEF's u10 `tanh` slightly exceeds
  1 ulp (1.0005).
