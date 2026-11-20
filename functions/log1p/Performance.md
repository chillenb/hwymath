# log1p (double precision): performance

Accuracy and throughput of the log1p implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512 and AArch64 hardware) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `log1pl` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts results that are not correctly rounded), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference. "special": results where the
  reference is NaN, infinite or zero and the result is not exactly that (Highway's contrib functions do not
  handle x outside their domain).
- **Throughput**: ns/value for 4096 random doubles in [-0.9, 10] (cache-resident), best of 7 runs of 2000
  repetitions (`make log1p-hwy`). Intel's libsvml, SLEEF and libmvec are called per vector.
- **Toolchain**: conda-forge GCC 15.3 (pixi), glibc 2.28 sysroot (old scalar libm symbol versions).
- The machine was shared with other jobs during the runs (load ~3 on 6 cores): timings vary by ~10 %.

## Summary

The Highway version is bit-identical to Intel's AVX-512 code on every FMA target, so its accuracy is that of
SVML's AVX-512 algorithm on every machine.

## Intel Core i5-9500 (AVX2), local

20,031,535 inputs. All SVML-algorithm rows are bit-identical (AVX2 checked against the C/AVX-512 translation
under Intel SDE).

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML | 4 | 0.6064 | 0.1793 % | 0 | 0 | 2.95 |
| Highway contrib `Log1p`, AVX2 |  |  | 2.4268 | 10.6941 % | 311826 | 1251731 | 1.69 |
| Highway SSE4 / SSSE3 / SSE2 | `std::log1p` | 2 | 0.8091 | 2.3817 % | 0 | 0 | 10.14 |
| Intel libsvml HA AVX2 (`__svml_log1p4_ha_l9`) |  | 4 | 0.7407 | 0.2783 % | 0 | 0 | 2.75 |
| SLEEF u10 AVX2 |  | 4 | 1.0000 | 0.8001 % | 0 | 0 | 3.69 |
| glibc `log1p` (scalar) |  | 1 | 0.8091 | 2.3817 % | 0 | 0 | 10.17 |

## AVX-512 (Intel SDE `-spr`, accuracy only)

2,031,535 inputs (`make log1p-avx512-test`); SDE timings are meaningless.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C SVML translation AVX-512 = Highway AVX3/AVX3_DL/AVX3_SPR | 0.6172 | 0.1821 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_log1p8_ha_z0`) | 0.6172 | 0.1821 % | 0 |
| SLEEF u10 AVX-512 | 1.0000 | 0.8101 % | 0 |

Intel's libsvml binary `__svml_log1p8_ha_z0` gives the same results as this version except for the bits of NaN results (x < -1): 0 differences in non-NaN results on 531,535 inputs under SDE.
