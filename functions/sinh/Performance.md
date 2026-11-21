# sinh (double precision): performance

Accuracy and throughput of the sinh implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result against `sinhl` (`long double`), over the
  inputs of `tests/inputs.h` (special values, ±50 doubles around every threshold, table boundaries,
  uniform values in [−720, 720] and [−25, 25], log-uniform magnitudes 2^−60 ... 2^4, random bit patterns);
  statistics over the inputs with a finite, nonzero reference, and exact agreement required for
  NaN/inf/zero results ("special" column: mismatches). "not CR": error > 0.5 ulp. The x86-64 `long double`
  has a 64-bit mantissa, which misjudges some results lying almost exactly halfway between two doubles, so
  the x86-64 "not CR" percentages are slightly optimistic; the AArch64 run (binary128 reference) shows the
  difference.
- **Throughput** (`make sinh-hwy`): nanoseconds per value for 4096 random doubles uniform in [−600, 600]
  (cache-resident, no callout lanes), best of 7 runs of 2000 repetitions, on a Core i5-9500 (Coffee Lake,
  AVX2, no AVX-512) shared with other jobs: expect a few percent of noise. Intel's libsvml and SLEEF are
  called per vector of 4 lanes through the same array loop. SDE timings are meaningless and not given;
  AVX-512 timings will come from CI.
- **Toolchain**: conda-forge GCC 15.3, Highway 1.4.0, SLEEF 3.9.0, Intel libsvml from
  `intel-cmplr-lib-rt` 2026.2.0 (comparison only). glibc's scalar `sinh` binds to the glibc 2.28 sysroot's symbol.

## Core i5-9500 (Coffee Lake, AVX2), local

`make sinh-hwy`: 20,007,045 inputs, x86-64 `long double` reference.

| implementation | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm) | 4 | **0.8696** | **0.5449 %** | 0 | **1.93** |
| Highway SSE4 / SSSE3 / SSE2 (`std::sinh`) | 2 | 1.8564 | 9.9407 % | 43,166 | 11.6–11.8 |
| Intel libsvml HA AVX2 (`__svml_sinh4_ha_l9`) | 4 | 0.8696 | 0.5449 % | 0 | 3.16 |
| SLEEF u10 AVX2 (`Sleef_sinhd4_u10avx2`) | 4 | inf* | 0.0947 % | 5,165* | 4.50 |
| glibc `sinh` (scalar) | 1 | 1.8564 | 9.9407 % | 43,166 | 11.88 |

All Highway targets with FMA give bit-identical results, so the accuracy figures of the Highway
AVX2 row hold for every such target (AVX-512, NEON, SVE).

## AVX-512 (Intel SDE, `-spr`), accuracy only

`make sinh-avx512-test`: 2,007,045 inputs.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C translation AVX-512 = Highway AVX3_SPR/AVX3_DL/AVX3/AVX2 (bit-identical) | 0.8359 | 0.5437 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_sinh8_ha_z0`) | 0.8359 | 0.5437 % | 0 |
| SLEEF u10 AVX-512 | inf* | 0.1081 % | 688* |

## AArch64 (QEMU), binary128 reference

`common/tools/xcheck.sh` (as in `make sinh-xcheck`, with `NR` = 300,000): 300,000 random inputs plus the special
and threshold inputs: NEON, SVE and SVE2 at 128 to 2048 bits are bit-identical to x86-64; against the
binary128 `long double` reference the error is max 0.8311 ulp, 0.5872 % not correctly rounded (the same
results measured with the x86-64 reference: 0.8306 ulp, 0.5715 %).

## Observations

- Highway's AVX2 version of Intel's algorithm is about 1.6x as fast as Intel's own AVX2 HA sinh in libsvml,
  which has the same error statistics on these inputs (probably the same algorithm). Its 4-lane table
  lookups by rows (README.md) matter: with the permute/blend tree it ran at 2.41–2.45 ns/value.
- Max error 0.87 ulp, below the 1 ulp HA bound; 0.54 % of the results are not correctly rounded.
- glibc's scalar `sinh` reaches 1.86 ulp (43,166 results above 1 ulp).
- \* SLEEF's u10 `sinh` returns inf for 709.78 < |x| < 710.476, where the result is still finite (its
  documented domain ends at 709); elsewhere it is the most often correctly rounded.
