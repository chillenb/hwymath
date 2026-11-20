# asinh (double precision): performance

Accuracy and throughput of the asinh implementations in this directory and of the libraries they are compared
with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `asinhl` (`long double`), over the
  test inputs of `tests/inputs.h` (20,000,000 random values plus special values and thresholds, see
  [README.md](README.md)); statistics over the inputs with a finite, nonzero result. "not CR": results that
  are not correctly rounded (error > 0.5 ulp). "special": results that differ from the reference where it is
  ±0, ±inf or NaN. The x86-64 `long double` has a 64-bit mantissa, so the "not CR" percentages are slightly
  optimistic.
- **Throughput**: nanoseconds per value for 4096 random doubles uniform in [−20, 20] (cache-resident), best of
  7 runs of 1000 repetitions (`make asinh-hwy`). Array functions process the whole array per call; Intel's
  libsvml and SLEEF are called once per vector through a small loop.
- **Toolchain**: conda-forge GCC 15.3 via pixi (glibc 2.28 sysroot: the scalar glibc baseline may bind to old
  symbol versions). The machine is shared with other jobs, so expect some variation in the timings (the
  accuracy figures are exact).

Intel SDE timings are meaningless; the AVX-512 numbers come from CI (added at integration). Every Highway
target with FMA is bit-identical to Intel's AVX-512 code, so the results are the same on every machine:
under Intel SDE the AVX3, AVX3_DL and AVX3_SPR targets and the C/AVX-512 translation give exactly the AVX2
results below, and so does Intel's libsvml `__svml_asinh8_ha_z0` (apart from NaN payloads). On AArch64 the
error statistics differ slightly only because the `long double` reference is more precise (binary128).

## Core i5-9500 (Coffee Lake, AVX2), local

20,000,000 random inputs plus specials (`./build/asinh/hwy_test`), `long double` with 64-bit mantissa.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML (bit-identical to `__svml_asinh8_ha`) | 4 | 0.7383 | 0.4969 % | 0 | 0 | 7.96 |
| Highway SSE4 / SSSE3 / SSE2 | `std::asinh` | 2 | 1.6899 | 12.5794 % | 24,874 | 0 | 11.4 |
| Intel SVML HA AVX2 (`__svml_asinh4_ha_l9`, a different implementation) | | 4 | 0.7134 | 0.1258 % | 0 | 0 | 6.70 |
| SLEEF u10 AVX2 | | 4 | inf | 15.7026 % | 3,124,475 | 0 | 11.69 |
| glibc `asinh` (scalar) | | 1 | 1.6899 | 12.5794 % | 24,874 | 0 | 11.51 |

Largest error of the SVML algorithm: at x = −0x1.ff73a21a46e1p−8 (≈ −0.0078, small-input polynomial).

SLEEF's u10 `asinh` returns ±inf for |x| ≳ 1.34·10¹⁵⁴ (x² overflows), hence its "inf" maximum.

## Cost of the exact emulations (AVX2, i5-9500)

The same Highway code with the emulated instructions replaced (only to measure; these variants are not
bit-identical and not part of the build):

| variant | ns/value |
|---|---|
| exact `vrsqrt14pd` and `vrcp14pd` emulations (the shipped version) | 8.1 |
| `vrsqrt14pd` → 1/√x (`Sqrt` + `Div`) | 5.9 |
| both → 1/√x and 1/x | 4.5 |
| neither (wrong results, lower bound for the rest of the algorithm) | 3.4 |

So the exact emulations cost about 3.5 ns/value here, mostly latency: the coefficients of the linear
approximation are computed with a division (rcp14) and a square root plus a division (rsqrt14), which sit
on the critical path. Without evaluating two vectors per iteration and GCC's pre-allocation scheduling
(`-fschedule-insns -fsched-pressure`, set in `function.mk`) the shipped version takes 11.5 ns/value.
On CPUs with faster dividers (Ice Lake and later, Zen) the gap should be smaller.
