# log (double precision): performance

Accuracy and throughput of the log implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512 and AArch64 hardware) are added
at integration.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `logl` (`long double`, 64-bit
  mantissa on x86-64, which slightly undercounts results that are not correctly rounded), over the inputs of
  `tests/inputs.h`; statistics over the inputs with a finite, nonzero reference. "special": results where the
  reference is NaN, infinite or zero and the result is not exactly that (Highway's contrib functions do not
  handle x outside their domain).
- **Throughput**: ns/value for 4096 random doubles in [0.001, 1000] (cache-resident), best of 7 runs of 2000
  repetitions (`make log-hwy`). Intel's libsvml, SLEEF and libmvec are called per vector.
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
| Highway AVX2 | SVML | 4 | 0.6118 | 0.7784 % | 0 | 0 | 2.59 |
| Highway contrib `Log`, AVX2 |  |  | 0.8745 | 2.5441 % | 0 | 2502787 | 1.33 |
| Highway SSE4 / SSSE3 / SSE2 | `std::log` | 2 | 0.5166 | 0.0236 % | 0 | 0 | 4.21 |
| Intel libsvml HA AVX2 (`__svml_log4_ha_l9`) |  | 4 | 0.7437 | 0.3156 % | 0 | 0 | 1.80 |
| SLEEF u10 AVX2 |  | 4 | 0.7510 | 0.8422 % | 0 | 0 | 3.33 |
| glibc libmvec AVX2 |  | 4 | 1.5850 | 20.3823 % | 68326 | 0 | 2.68 |
| glibc `log` (scalar) |  | 1 | 0.5166 | 0.0236 % | 0 | 0 | 4.42 |

## AVX-512 (Intel SDE `-spr`, accuracy only)

2,029,320 inputs (`make log-avx512-test`); SDE timings are meaningless.

| implementation | max ulp | not CR | > 1 ulp |
|---|---|---|---|
| C SVML translation AVX-512 = Highway AVX3/AVX3_DL/AVX3_SPR | 0.6064 | 0.7690 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_log8_ha_z0`) | 0.5542 | 0.3646 % | 0 |
| SLEEF u10 AVX-512 | 0.7261 | 0.8435 % | 0 |
| glibc libmvec AVX-512 | 1.5884 | 20.3965 % | 6595 |

Intel's libsvml binary (conda-forge `intel-cmplr-lib-rt`) is a newer log than the studied source: 4,623 of the 500,000+ SDE results differ.
