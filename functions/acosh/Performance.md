# acosh (double precision): performance

Accuracy and throughput of the acosh implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `acoshl` (`long double`), over the
  test inputs of `tests/inputs.h` (20,000,000 random values plus special values and thresholds, see
  [README.md](README.md)); statistics over the inputs with a finite, nonzero result. "not CR": results that
  are not correctly rounded (error > 0.5 ulp). "special": results that differ from the reference where it is
  ±0, ±inf or NaN. The x86-64 `long double` has a 64-bit mantissa, so the "not CR" percentages are slightly
  optimistic.
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [1, 20] (cache-resident), best of
  7 runs of 1000 repetitions (`make acosh-hwy`). Array functions process the whole array per call; Intel's
  libsvml and SLEEF are called once per vector through a small loop.
- **Toolchain**: conda-forge GCC 15.3 via pixi (glibc 2.28 sysroot: the scalar glibc baseline may bind to old
  symbol versions). The machine is shared with other jobs, so expect some variation in the timings (the
  accuracy figures are exact).

Intel SDE timings are meaningless; the AVX-512 numbers come from CI (added at integration). Every Highway
target with FMA is bit-identical to Intel's AVX-512 code, so the results are the same on every machine:
under Intel SDE the AVX3, AVX3_DL and AVX3_SPR targets and the C/AVX-512 translation give exactly the AVX2
results below, and so does Intel's libsvml `__svml_acosh8_ha_z0` (apart from NaN payloads). On AArch64 the
error statistics differ slightly only because the `long double` reference is more precise (binary128).

## Core i5-9500 (Coffee Lake, AVX2), local

20,000,000 random inputs plus specials (`./build/acosh/hwy_test`), `long double` with 64-bit mantissa.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML (bit-identical to `__svml_acosh8_ha`) | 4 | 0.7349 | 0.2593 % | 0 | 0 | 8.79 |
| Highway SSE4 / SSSE3 / SSE2 | `std::acosh` | 2 | 2.1899 | 17.7307 % | 67,527 | 0 | 9.4 |
| Intel SVML HA AVX2 (`__svml_acosh4_ha_l9`, a different implementation) | | 4 | 0.8999 | 2.4882 % | 0 | 0 | 6.10 |
| SLEEF u10 AVX2 | | 4 | inf | 17.2538 % | 1,874,233 | 0 | 9.60 |
| glibc `acosh` (scalar) | | 1 | 2.1899 | 17.7307 % | 67,527 | 0 | 9.38 |

Largest error of the SVML algorithm: at x = 0x1.0001f4e6a6d75p+0 (≈ 1.00003, small-input polynomial).

SLEEF's u10 `acosh` returns +inf for x ≳ 1.34·10¹⁵⁴ (x² overflows), hence its "inf" maximum.

The exact `vrsqrt14pd`/`vrcp14pd` emulations cost about 3.5 ns/value on this machine (measured for asinh,
which has the same structure: see [../asinh/Performance.md](../asinh/Performance.md)).
