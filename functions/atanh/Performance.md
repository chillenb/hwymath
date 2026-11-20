# atanh (double precision): performance

Accuracy and throughput of the atanh implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `atanhl` (`long double`), over the
  test inputs of `tests/inputs.h` (20,000,000 random values plus special values and thresholds, see
  [README.md](README.md)); statistics over the inputs with a finite, nonzero result. "not CR": results that
  are not correctly rounded (error > 0.5 ulp). "special": results that differ from the reference where it is
  ±0, ±inf or NaN. The x86-64 `long double` has a 64-bit mantissa, so the "not CR" percentages are slightly
  optimistic.
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [−1, 1] (cache-resident), best of
  7 runs of 1000 repetitions (`make atanh-hwy`). Array functions process the whole array per call; Intel's
  libsvml and SLEEF are called once per vector through a small loop.
- **Toolchain**: conda-forge GCC 15.3 via pixi (glibc 2.28 sysroot: the scalar glibc baseline may bind to old
  symbol versions). The machine is shared with other jobs, so expect some variation in the timings (the
  accuracy figures are exact).

Intel SDE timings are meaningless; the AVX-512 numbers come from CI (added at integration). Every Highway
target with FMA is bit-identical to Intel's AVX-512 code, so the results are the same on every machine:
under Intel SDE the AVX3, AVX3_DL and AVX3_SPR targets and the C/AVX-512 translation give exactly the AVX2
results below, and so does Intel's libsvml `__svml_atanh8_ha_z0` (apart from NaN payloads). On AArch64 the
error statistics differ slightly only because the `long double` reference is more precise (binary128).

## Core i5-9500 (Coffee Lake, AVX2), local

20,000,000 random inputs plus specials (`./build/atanh/hwy_test`), `long double` with 64-bit mantissa.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML (bit-identical to `__svml_atanh8_ha`) | 4 | 0.7178 | 0.0117 % | 0 | 0 | 5.99 |
| Highway SSE4 / SSSE3 / SSE2 | `std::atanh` | 2 | 1.6118 | 7.0350 % | 31,599 | 0 | 18.0 |
| Intel SVML HA AVX2 (`__svml_atanh4_ha_l9`, a different implementation) | | 4 | 0.7129 | 0.0324 % | 0 | 0 | 5.66 |
| SLEEF u10 AVX2 | | 4 | 1.0000 | 0.0517 % | 0 | 0 | 5.81 |
| glibc `atanh` (scalar) | | 1 | 1.6118 | 7.0350 % | 31,599 | 0 | 17.68 |

Largest error of the SVML algorithm: at x = 0x1.0826222b1e24ep−6 (≈ 0.0161).

The two exact `vrcp14pd` emulations (one division each) cost about 1.3 ns/value each on this machine
(estimated from the asinh measurements in [../asinh/Performance.md](../asinh/Performance.md)).
