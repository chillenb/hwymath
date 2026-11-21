# cosh (double precision): performance

Accuracy and throughput of the cosh implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result against `coshl` (`long double`), over the
  inputs of `tests/inputs.h` (special values, ±50 doubles around every threshold, table boundaries,
  uniform values in [−720, 720] and [−25, 25], log-uniform magnitudes 2^−60 ... 2^4, random bit patterns);
  statistics over the inputs with a finite, nonzero reference, and exact agreement required for
  NaN/inf/zero results ("special" column: mismatches). "not CR": error > 0.5 ulp. The x86-64 `long double`
  has a 64-bit mantissa, which misjudges some results lying almost exactly halfway between two doubles, so
  the x86-64 "not CR" percentages are slightly optimistic; the AArch64 run (binary128 reference) shows the
  difference.
- **Throughput** (`make cosh-hwy`): nanoseconds per value for 4096 random doubles uniform in [−600, 600]
  (cache-resident, no callout lanes), best of 7 runs of 2000 repetitions, on a Core i5-9500 (Coffee Lake,
  AVX2, no AVX-512) shared with other jobs: expect a few percent of noise. Intel's libsvml and SLEEF are
  called per vector of 4 lanes through the same array loop. SDE timings are meaningless and not given;
  AVX-512 timings will come from CI.
- **Toolchain**: conda-forge GCC 15.3, Highway 1.4.0, SLEEF 3.9.0, Intel libsvml from
  `intel-cmplr-lib-rt` 2026.2.0 (comparison only). glibc's scalar `cosh` binds to the glibc 2.28 sysroot's symbol.

## Core i5-9500 (Coffee Lake, AVX2), local

`make cosh-hwy`: 20,006,641 inputs, x86-64 `long double` reference.

| implementation | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm) | 4 | **0.7510** | **3.2851 %** | 0 | **1.31** |
| Highway SSE4 / SSSE3 / SSE2 (`std::cosh`) | 2 | 1.8135 | 6.1841 % | 1,081 | 7.1–7.2 |
| Intel libsvml HA AVX2 (`__svml_cosh4_ha_l9`) | 4 | 0.7500 | 2.7612 % | 0 | 2.47 |
| SLEEF u10 AVX2 (`Sleef_coshd4_u10avx2`) | 4 | inf* | 0.0348 % | 5,165* | 4.45 |
| glibc `cosh` (scalar) | 1 | 1.8135 | 6.1841 % | 1,081 | 7.23 |

All Highway targets with FMA give bit-identical results, so the accuracy figures of the Highway
AVX2 row hold for every such target (AVX-512, NEON, SVE).

## AVX-512 (Intel SDE, `-spr`), accuracy only

`make cosh-avx512-test`: 2,006,641 inputs.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C translation AVX-512 = Highway AVX3_SPR/AVX3_DL/AVX3/AVX2 (bit-identical) | 0.7505 | 3.3010 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_cosh8_ha_z0`) | 0.7505 | 3.3010 % | 0 |
| SLEEF u10 AVX-512 | inf* | 0.0451 % | 688* |

## AArch64 (QEMU), binary128 reference

`common/tools/xcheck.sh` (as in `make cosh-xcheck`, with `NR` = 300,000): 300,000 random inputs plus the special
and threshold inputs: NEON, SVE and SVE2 at 128 to 2048 bits are bit-identical to x86-64; against the
binary128 `long double` reference the error is max 0.7502 ulp, 3.3587 % not correctly rounded (the same
results measured with the x86-64 reference: 0.7500 ulp, 3.3416 %).

## Observations

- Highway's AVX2 version is about 1.9x as fast as Intel's AVX2 HA cosh in libsvml, which uses a
  different (slightly more often correctly rounded) algorithm on AVX2; Intel's AVX-512 libsvml cosh has
  exactly the statistics of the open-source algorithm ported here. With the permute/blend-tree lookups the
  Highway version ran at 1.72–1.74 ns/value.
- Max error 0.75 ulp, below the 1 ulp HA bound; 3.3 % of the results are not correctly rounded.
- \* SLEEF's u10 `cosh` returns inf for 709.78 < |x| < 710.476 (documented domain |x| ≤ 709).
