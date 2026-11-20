# cbrt (double precision): performance

Accuracy and throughput of the cbrt implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md).

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `cbrtl` (`long double`), over the
  test inputs of `tests/inputs.h`: special values, powers of two over the whole exponent range, perfect
  cubes, the algorithm's thresholds with neighbouring doubles, random bit patterns, random mantissas with
  small exponents, and random values in [−1000, 1000]. The statistics cover the inputs with a finite,
  nonzero result; 0, inf and NaN must match exactly ("special" column). "not CR" counts results that are
  not correctly rounded (error > 0.5 ulp). The x86-64 `long double` has a 64-bit mantissa, which slightly
  undercounts not-correctly-rounded results.
- **Throughput**: nanoseconds per value for 4096 random doubles in [−1000, 1000] (cache-resident), best of
  7 runs of 2000 repetitions, array functions (Intel's and SLEEF's vector functions are called once per
  vector).
- **Toolchain**: conda-forge GCC 15.3 via pixi, glibc 2.28 sysroot (scalar `cbrt` binds to that
  version).

Program: `make cbrt-hwy` (`tests/hwy_test.cc`).

## Core i5-9500 (Coffee Lake, AVX2), local

20,019,282 inputs.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | ns/value |
|---|---|---|---|---|---|---|
| **Highway AVX2** | SVML, bit-identical | 4 | **0.5547** | 0.5994 % | 0 | **2.23** |
| Highway SSE4 / SSSE3 / SSE2 | `std::cbrt` (no FMA) | 2 | 3.4126 | 52.3089 % | 4,043,782 | 10.7 |
| Intel libsvml HA AVX2 (`__svml_cbrt4_ha_l9`) | | 4 | 0.5156 | 0.1532 % | 0 | 2.27 |
| SLEEF u10 AVX2 | | 4 | 0.5410 | 0.0308 % | 0 | 6.63 |
| C translation of the SVML scalar routine | | 1 | 0.5015 | 0.0038 % | 0 | 36.52 |
| glibc `cbrt` (scalar) | | 1 | 3.4126 | 52.3089 % | 4,043,782 | 10.76 |

Highway AVX2 is as fast as Intel's own AVX2 HA routine (a different algorithm), 3x faster than SLEEF u10
and 5x faster than glibc's scalar `cbrt`, and it is far more accurate than glibc's.

## AVX-512 (Intel SDE `-spr` on the same machine; accuracy only)

2,019,282 inputs. SDE timings are meaningless; CI numbers are added at integration.

| implementation | lanes | max ulp | not CR | > 1 ulp |
|---|---|---|---|---|
| Highway AVX3_SPR / AVX3_DL / AVX3 / AVX2 (bit-identical to each other and to the C translation) | 8 / 4 | 0.5469 | 0.6057 % | 0 |
| C translation of `__svml_cbrt8_ha` | 8 | 0.5469 | 0.6057 % | 0 |
| Intel libsvml HA AVX-512 (`__svml_cbrt8_ha_z0`) | 8 | 0.5469 | 0.6057 % | 0 |
| Intel libsvml HA AVX2 (`__svml_cbrt4_ha_l9`) | 4 | 0.5137 | 0.1545 % | 0 |
| SLEEF u10 AVX2 / AVX-512 | 4 / 8 | 0.5400 | 0.0934 % | 0 |
| glibc `cbrt` (scalar) | 1 | 3.2842 | 52.2720 % | 406,432 |

libsvml's `__svml_cbrt8_ha_z0` gives the same statistics as the translated `__svml_cbrt8_ha`. It is
presumably the same code.
