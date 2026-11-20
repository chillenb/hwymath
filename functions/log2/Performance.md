# log2 (double precision): performance

Accuracy and throughput of the log2 implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512 and AArch64 hardware) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `log2l` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts results that are not correctly rounded), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference. "special": results where the
  reference is NaN, infinite or zero and the result is not exactly that (Highway's contrib functions do not
  handle x outside their domain).
- **Throughput**: ns/value for 4096 random doubles in [0.001, 1000] (cache-resident), best of 7 runs of 2000
  repetitions (`make log2-hwy`). Intel's libsvml, SLEEF and libmvec are called per vector.
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
| Highway AVX2 | SVML | 4 | 0.5195 | 0.0523 % | 0 | 0 | 2.53 |
| Highway contrib `Log2`, AVX2 |  |  | 1.6538 | 28.7032 % | 422764 | 2502787 | 1.39 |
| Highway SSE4 / SSSE3 / SSE2 | `std::log2` | 2 | 0.5444 | 0.0687 % | 0 | 0 | 4.56 |
| Intel libsvml HA AVX2 (`__svml_log24_ha_l9`) |  | 4 | 0.5024 | 0.0070 % | 0 | 0 | 2.28 |
| SLEEF u10 AVX2 |  | 4 | 0.8262 | 2.5023 % | 0 | 0 | 3.47 |
| glibc `log2` (scalar) |  | 1 | 0.5444 | 0.0687 % | 0 | 0 | 4.46 |

## AVX-512 (Intel SDE `-spr`, accuracy only)

2,029,320 inputs (`make log2-avx512-test`); SDE timings are meaningless.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C SVML translation AVX-512 = Highway AVX3/AVX3_DL/AVX3_SPR | 0.5137 | 0.0501 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_log28_ha_z0`) | 0.5137 | 0.0501 % | 0 |
| SLEEF u10 AVX-512 | 0.8096 | 2.4927 % | 0 |

Intel's libsvml binary `__svml_log28_ha_z0` gives exactly the same results as this version (0 differences on 529,320 inputs under SDE).
