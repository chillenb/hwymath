# acos (double precision): performance

Accuracy and throughput of the acos implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (other machines) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `acosl` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts not-correctly-rounded results), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference; "special" counts
  NaN/inf/signed-zero mismatches where the reference is not finite or zero (none anywhere).
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [−1, 1] (cache-resident), best
  of 7 runs of 2000 repetitions. Highway versions are array functions; libsvml and SLEEF are called once
  per 4- or 8-lane vector.
- Program: `make acos-hwy` (`tests/hwy_test.cc`). Toolchain: conda-forge GCC via pixi (glibc 2.28 sysroot).

## Core i5-9500 (AVX2), local, 20,004,655 inputs

| implementation | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm, portable rsqrt) | 0.8174 | 6.7724 % | 0 | **3.06** |
| Highway SSE4/SSSE3/SSE2 (`std::acos`) | 0.5210 | 0.0370 % | 0 | 12.7–12.8 |
| Intel SVML HA AVX2 (libsvml `__svml_acos4_ha_l9`) | 0.8174 | 6.7690 % | 0 | 3.00 |
| SLEEF u10 AVX2 | 0.7778 | 1.7303 % | 0 | 4.04 |
| glibc `acos` (scalar) | 0.5210 | 0.0370 % | 0 | 12.84 |

The machine is shared with other jobs: expect some run-to-run variation in the timings.

## AVX-512 targets under Intel SDE (`-spr`), 2,004,655 inputs (accuracy only)

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (portable rsqrt; identical results) | 0.8174 | 6.7737 % | 0 |
| Highway AVX3* native `vrsqrt14pd` (= Intel's results) | 0.8174 | 6.7588 % | 0 |
| C translation `svml_acos8_ha` (= `__svml_acos8_ha`) | 0.8174 | 6.7588 % | 0 |
| C translation, portable rsqrt | 0.8174 | 6.7737 % | 0 |
| Intel SVML HA AVX-512 (libsvml `__svml_acos8_ha_z0`) | 0.8174 | 6.7679 % | 0 |
| Intel SVML HA AVX2 (libsvml `__svml_acos4_ha_l9`) | 0.8115 | 6.7705 % | 0 |
| SLEEF u10 AVX-512 / AVX2 | 0.7383 | 1.7231 % | 0 |
| glibc `acos` (scalar) | 0.5186 | 0.0403 % | 0 |

The portable rsqrt approximation changes 483 of the 2,004,655 results relative to Intel's (0.024 %).
SDE timings are meaningless and not shown.
