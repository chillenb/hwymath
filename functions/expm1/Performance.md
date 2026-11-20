# expm1 (double precision): performance

Accuracy and throughput of the expm1 implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512, AArch64)
are added at integration. Timings under Intel SDE are meaningless and are not reported.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `expm1l` (`long double`), over
  the test inputs of `tests/inputs.h` (20,014,924 inputs: special values, every threshold with
  neighbouring doubles, table boundaries, uniform in [−750, 750] and [−1, 1], random bit patterns).
  The statistics cover the inputs with a finite, nonzero result. "not CR" means not correctly rounded
  (error > 0.5 ulp). "special" counts NaN/inf/zero mismatches. The x86-64 `long double` has a 64-bit
  mantissa, so the x86-64 "not CR" percentages are slightly optimistic. With AArch64's binary128
  (QEMU, `make expm1-xcheck`, 2M inputs), the SVML algorithm has max 0.5204 ulp and 0.0132 % not CR.
- **Throughput**: ns/value for 4096 random doubles in [−700, 700] (cache-resident, no callouts),
  best of 7 runs of 2000 repetitions, array functions (`make expm1-hwy`).
- **Toolchain**: conda-forge GCC via pixi, Highway 1.4, glibc 2.28 sysroot.

## Core i5-9500 (AVX2), local

Accuracy: `hwy_test --no-bench` (20M random inputs). Timings: `hwy_test 2000000` on a shared machine,
so expect some noise.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| **Highway AVX2** | SVML (bit-identical) | 4 | **0.5264** | 0.0083 % | 0 | 0 | **1.95** |
| Highway SSE4/SSSE3/SSE2 | `std::expm1` | 2 | 0.8799 | 3.9407 % | 0 | 0 | 7.87 |
| Highway contrib `Expm1`, AVX2 (x ≤ 706) | | | 4.0127 | 10.1556 % | 268972 | 1 | 1.68 |
| Intel SVML HA AVX2 (`__svml_expm14_ha_l9`) | | | 0.6060 | 0.0210 % | 0 | 0 | 3.07 |
| SLEEF u10 AVX2 | | | 1.0000 | 0.0394 % | 0 | 0 | 5.13 |
| glibc `expm1` (scalar) | | 1 | 0.8799 | 3.9407 % | 0 | 0 | 7.81 |

## AVX-512 (accuracy only, Intel SDE `-spr`, 2,014,924 inputs)

All Highway FMA targets (AVX3_SPR, AVX3_DL, AVX3, AVX2) give results bit-identical to the C/AVX-512
translation, and so to `__svml_expm18_ha`: max 0.5200 ulp, 0.0095 % not CR, none > 1 ulp. Intel's
libsvml `__svml_expm18_ha_z0` measures the same (0.5200 ulp, 0.0095 %). SLEEF u10 AVX-512 has max
1.0000 ulp.

## Error by path (Highway AVX2, 20M inputs)

| lanes | max ulp | not CR |
|---|---|---|
| vector path, \|x\| < 1 | 0.5264 | 0.0143 % |
| vector path, \|x\| ≥ 1 (x < 707.7) | 0.5010 | 0.0001 % |
| callout (707.7 ≤ x < 709.78) | 0.5054 | 0.0107 % |
