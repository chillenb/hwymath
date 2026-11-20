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
| **Highway AVX2** (SVML algorithm, portable reciprocal) | 0.7920 | 2.8815 % | 0 | **4.06** |
| Highway SSE4/SSSE3/SSE2 (`std::atan2`) | 0.5225 | 0.0332 % | 0 | 23.0–23.8 |
| Intel SVML HA AVX2 (libsvml `__svml_atan24_ha_l9`) | 0.7920 | 2.8569 % | 0 | 5.02 |
| SLEEF u10 AVX2 | 0.9932 | 1.2660 % | 0 | 6.08 |
| glibc `atan2` (scalar) | 0.5225 | 0.0332 % | 0 | 23.00 |

The machine is shared with other jobs: expect some run-to-run variation in the timings.

## AVX-512 targets under Intel SDE (`-spr`), 2,003,640 input pairs (accuracy only)

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (portable reciprocal; identical results) | 0.7749 | 2.8922 % | 0 |
| Highway AVX3* native `vrcp14pd` (= Intel's results) | 0.7749 | 2.8921 % | 0 |
| C translation `svml_atan28_ha` (= `__svml_atan28_ha`) | 0.7749 | 2.8921 % | 0 |
| C translation, portable reciprocal | 0.7749 | 2.8922 % | 0 |
| Intel SVML HA AVX-512 (libsvml `__svml_atan28_ha_z0`) | 0.7749 | 2.8660 % | 0 |
| Intel SVML HA AVX2 (libsvml `__svml_atan24_ha_l9`) | 0.7749 | 2.8656 % | 0 |
| SLEEF u10 AVX-512 / AVX2 | 0.9395 | 1.2674 % | 0 |
| glibc `atan2` (scalar) | 0.5210 | 0.0347 % | 0 |

The portable reciprocal changes 6 of the 2,003,640 results relative to Intel's (0.0003 %). The
libsvml binaries are a slightly different build than the open-source assembly (their not-CR rate is
0.026 points lower). SDE timings are meaningless and not shown.
