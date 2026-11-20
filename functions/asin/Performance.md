# asin (double precision): performance

Accuracy and throughput of the asin implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (other machines) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `asinl` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts not-correctly-rounded results), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference; "special" counts
  NaN/inf/signed-zero mismatches where the reference is not finite or zero (none anywhere).
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [−1, 1] (cache-resident), best
  of 7 runs of 2000 repetitions. Highway versions are array functions; libsvml and SLEEF are called once
  per 4- or 8-lane vector.
- Program: `make asin-hwy` (`tests/hwy_test.cc`). Toolchain: conda-forge GCC via pixi (glibc 2.28 sysroot).

## Core i5-9500 (AVX2), local, 20,004,655 inputs

| implementation | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm, portable rsqrt) | 0.8389 | 6.1196 % | 0 | **2.31** |
| Highway SSE4/SSSE3/SSE2 (`std::asin`) | 0.5146 | 0.1007 % | 0 | 11.9–12.4 |
| Intel SVML HA AVX2 (libsvml `__svml_asin4_ha_l9`) | 0.8389 | 6.1136 % | 0 | 2.36 |
| SLEEF u10 AVX2 | 0.8984 | 4.2965 % | 0 | 3.11 |
| glibc `asin` (scalar) | 0.5146 | 0.1007 % | 0 | 12.20 |

The machine is shared with other jobs: expect some run-to-run variation in the timings.

## AVX-512 targets under Intel SDE (`-spr`), 2,004,655 inputs (accuracy only)

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (portable rsqrt; identical results) | 0.8291 | 6.1020 % | 0 |
| Highway AVX3* native `vrsqrt14pd` (= Intel's results) | 0.8291 | 6.0793 % | 0 |
| C translation `svml_asin8_ha` (= `__svml_asin8_ha`) | 0.8291 | 6.0793 % | 0 |
| C translation, portable rsqrt | 0.8291 | 6.1020 % | 0 |
| Intel SVML HA AVX-512 (libsvml `__svml_asin8_ha_z0`) | 0.8291 | 6.0919 % | 0 |
| Intel SVML HA AVX2 (libsvml `__svml_asin4_ha_l9`) | 0.8223 | 6.0927 % | 0 |
| SLEEF u10 AVX-512 / AVX2 | 0.8618 | 4.2940 % | 0 |
| glibc `asin` (scalar) | 0.5103 | 0.0987 % | 0 |

The portable rsqrt approximation changes 620 of the 2,004,655 results relative to Intel's (0.031 %).
The libsvml binaries are a slightly different build than the open-source assembly. SDE timings are
meaningless and not shown.
