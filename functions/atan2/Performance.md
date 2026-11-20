# atan2 (double precision): performance

Accuracy and throughput of the atan2 implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (other machines) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `atan2l` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts not-correctly-rounded results), over the input pairs of
  `tests/inputs.h`; statistics over the pairs with a finite, nonzero reference; "special" counts NaN/inf/
  signed-zero mismatches where the reference is not finite or zero (none anywhere).
- **Throughput**: nanoseconds per value for 4096 random pairs, y and x uniform in [-2, 2] (cache-resident),
  best of 7 runs of 1000 repetitions, array functions.
- Program: `make atan2-hwy` (`tests/hwy_test.cc`). Toolchain: conda-forge GCC via pixi (glibc 2.28 sysroot).

## Core i5-9500 (AVX2), local, 20,003,640 input pairs

| implementation | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|
| **Highway AVX2** (SVML algorithm, exact `vrcp14pd` emulation) | 0.7769 | 2.8849 % | 0 | **5.67** |
| Highway SSE4/SSSE3/SSE2 (`std::atan2`) | 0.5225 | 0.0344 % | 0 | 21.8–23.7 |
| Intel SVML HA AVX2 (libsvml `__svml_atan24_ha_l9`) | 0.7769 | 2.8601 % | 0 | 4.82 |
| SLEEF u10 AVX2 | 0.9932 | 1.2620 % | 0 | 5.76 |
| glibc `atan2` (scalar) | 0.5225 | 0.0344 % | 0 | 21.86 |

The machine is shared with other jobs: expect some run-to-run variation in the timings. The exact
`vrcp14pd` emulation (AVX2 only; AVX-512 targets use the instruction) is the price of bit-identity with
Intel: with the earlier, non-identical trunc16 substitute the AVX2 version took 4.09 ns/value.

## AVX-512 targets under Intel SDE (`-spr`), 2,003,640 input pairs (accuracy only)

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (bit-identical to Intel's code) | 0.7861 | 2.9002 % | 0 |
| C translation `svml_atan28_ha` (= `__svml_atan28_ha`) | 0.7861 | 2.9002 % | 0 |
| Intel SVML HA AVX-512 (libsvml `__svml_atan28_ha_z0`) | 0.7861 | 2.8734 % | 0 |
| Intel SVML HA AVX2 (libsvml `__svml_atan24_ha_l9`) | 0.7861 | 2.8730 % | 0 |
| SLEEF u10 AVX-512 / AVX2 | 0.9453 | 1.2719 % | 0 |
| glibc `atan2` (scalar) | 0.5210 | 0.0340 % | 0 |

The libsvml binaries are a slightly different build than the open-source assembly (their not-CR rate
is 0.027 points lower). SDE timings are meaningless and not shown.
