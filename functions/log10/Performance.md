# log10 (double precision): performance

Accuracy and throughput of the log10 implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512 and AArch64 hardware) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `log10l` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts results that are not correctly rounded), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference. "special": results where the
  reference is NaN, infinite or zero and the result is not exactly that (Highway's contrib functions do not
  handle x outside their domain).
- **Throughput**: ns/value for 4096 random doubles in [0.001, 1000] (cache-resident), best of 7 runs of 2000
  repetitions (`make log10-hwy`). Intel's libsvml, SLEEF and libmvec are called per vector.
- **Toolchain**: conda-forge GCC 15.3 (pixi), glibc 2.28 sysroot (old scalar libm symbol versions).
- The machine was shared with other jobs during the runs (load ~3 on 6 cores): timings vary by ~10 %.

## Summary

The Highway version is bit-identical to Intel's AVX-512 code on every FMA target, so its accuracy is that of
SVML's AVX-512 algorithm on every machine.

## Intel Core i5-9500 (AVX2), local

20,029,320 inputs. All SVML-algorithm rows are bit-identical (AVX2 checked against the C/AVX-512 translation
under Intel SDE).

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML | 4 | 0.6050 | 0.8371 % | 0 | 0 | 2.46 |
| Highway contrib `Log10`, AVX2 |  |  | 1.8857 | 29.3045 % | 481240 | 2502787 | 1.33 |
| Highway SSE4 / SSSE3 / SSE2 | `std::log10` | 2 | 1.5923 | 14.1384 % | 159537 | 0 | 7.96 |
| Intel libsvml HA AVX2 (`__svml_log104_ha_l9`) |  | 4 | 0.5107 | 0.0516 % | 0 | 0 | 2.42 |
| SLEEF u10 AVX2 |  | 4 | 0.7544 | 2.4826 % | 0 | 0 | 3.27 |
| glibc `log10` (scalar) |  | 1 | 1.5923 | 14.1384 % | 159537 | 0 | 8.25 |

## AVX-512 (Intel SDE `-spr`, accuracy only)

2,029,320 inputs (`make log10-avx512-test`); SDE timings are meaningless.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C SVML translation AVX-512 = Highway AVX3/AVX3_DL/AVX3_SPR | 0.6016 | 0.8378 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_log108_ha_z0`) | 0.6016 | 0.8378 % | 0 |
| SLEEF u10 AVX-512 | 0.7446 | 2.4716 % | 0 |

Intel's libsvml binary `__svml_log108_ha_z0` gives exactly the same results as this version (0 differences on 529,320 inputs under SDE).
