# atan (double precision): performance

Accuracy and throughput of the atan implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (other machines) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `atanl` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts not-correctly-rounded results), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference; "special" counts
  NaN/inf/signed-zero mismatches where the reference is not finite or zero (none anywhere).
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [−8, 8] (cache-resident), best
  of 7 runs of 2000 repetitions. Highway versions are array functions; libsvml and SLEEF are called once
  per 4- or 8-lane vector.
- Program: `make atan-hwy` (`tests/hwy_test.cc`). Toolchain: conda-forge GCC via pixi (glibc 2.28 sysroot).

## Core i5-9500 (AVX2), local, 20,004,201 inputs

| implementation | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm, exact `vrcp14pd` emulation) | 0.5640 | 0.0293 % | 0 | **6.29** |
| Highway SSE4/SSSE3/SSE2 (`std::atan`) | 0.5190 | 0.0349 % | 0 | 9.4–10.2 |
| Intel SVML HA AVX2 (libsvml `__svml_atan4_ha_l9`, a different algorithm) | 0.7324 | 3.4689 % | 0 | 4.10 |
| SLEEF u10 AVX2 | 0.8242 | 1.0018 % | 0 | 4.14 |
| glibc `atan` (scalar) | 0.5190 | 0.0349 % | 0 | 9.88 |

The machine is shared with other jobs: expect some run-to-run variation in the timings. The exact
`vrcp14pd` emulation (AVX2 only; AVX-512 targets use the instruction) is the price of bit-identity with
Intel: with the earlier, non-identical trunc16 substitute the AVX2 version took 4.51 ns/value.

## AVX-512 targets under Intel SDE (`-spr`), 2,004,201 inputs (accuracy only)

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (bit-identical to Intel's code) | 0.5645 | 0.0324 % | 0 |
| C translation `svml_atan8_ha` (= `__svml_atan8_ha`) | 0.5645 | 0.0324 % | 0 |
| Intel SVML HA AVX-512 (libsvml `__svml_atan8_ha_z0`) | 0.5645 | 0.0324 % | 0 |
| Intel SVML HA AVX2 (libsvml `__svml_atan4_ha_l9`) | 0.7227 | 3.4616 % | 0 |
| SLEEF u10 AVX-512 / AVX2 | 0.8135 | 1.0088 % | 0 |
| glibc `atan` (scalar) | 0.5151 | 0.0369 % | 0 |

SDE timings are meaningless and not shown.
