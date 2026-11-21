# Performance

Accuracy and speed are documented per function (`functions/<f>/Performance.md`, with the method and the
machines). This page summarizes the Highway versions as measured in CI.

Every Highway version except exp's is a port of Intel's own algorithm and returns bit-identical results to
Intel's AVX-512 code on every target with FMA; exp's Highway version is the more accurate "hp" algorithm
(a Highway port of Intel's exp algorithm is in progress).

## CI, ns per value (lower is better)

x86-64: Intel Xeon 6975P-C (Granite Rapids); AArch64: Neoverse-N2 (128-bit SVE2). GitHub Actions run
[44857079537](https://github.com/chillenb/hwymath/actions/runs/44857079537) (PR #10), shared machines: differences
below ~10 % are noise. Max error: worst case found over the ~20M test inputs, in ulps,
against a `long double` reference; identical on all FMA targets. Trig functions: |x| ≤ 100.

| function | max error (ulp) | Highway AVX-512 | Intel HA AVX-512 | Highway AVX2 | Intel HA AVX2 | SLEEF u10 AVX-512 | Highway N2 (best target) | SLEEF u10 AdvSIMD | glibc scalar (N2) |
|---|---|---|---|---|---|---|---|---|---|
| [exp](functions/exp/Performance.md) ¹ | 0.5107 | 0.55 | 0.39 | 0.90 | 0.84 | 0.67 | 3.79 | 4.59 | 7.02 |
| [exp2](functions/exp2/Performance.md) | 0.7520 | 0.53 | 0.37 | 1.60 | 0.83 | 0.63 | 4.84 | — | 7.05 |
| [expm1](functions/expm1/Performance.md) | 0.5264 | 0.76 | 0.76 | 2.02 | 2.38 | 3.00 | 6.01 | — | 8.01 |
| [log](functions/log/Performance.md) | 0.6118 | 0.73 | 0.64 | 2.63 | 1.37 | 1.11 | 7.73 | 9.66 | 4.15 |
| [log2](functions/log2/Performance.md) | 0.5195 | 0.76 | 0.70 | 2.69 | 1.48 | 1.17 | 7.77 | 10.29 | 4.70 |
| [log10](functions/log10/Performance.md) | 0.6050 | 0.75 | 0.70 | 2.69 | 1.53 | 1.10 | 7.86 | 10.40 | 7.63 |
| [log1p](functions/log1p/Performance.md) | 0.6064 | 1.01 | 0.83 | 3.03 | 1.96 | 1.50 | 9.31 | 9.47 | 11.28 |
| [pow](functions/pow/Performance.md) ² | 0.7505 | 4.04 | 1.46 | 10.84 | 4.31 | 4.44 | 30.90 | 26.26 | 13.23 |
| [cbrt](functions/cbrt/Performance.md) | 0.5547 | 1.00 | 0.78 | 2.13 | 1.85 | 2.39 | 8.25 | 15.99 | 13.94 |
| [sin](functions/sin/Performance.md) | 0.5200 | 0.73 | 0.78 | 2.10 | 1.93 | 1.97 | 4.49 | 10.80 | 13.24 |
| [cos](functions/cos/Performance.md) | 0.5195 | 0.70 | 0.70 | 1.99 | 1.93 | 2.20 | 4.18 | 12.16 | 14.57 |
| [tan](functions/tan/Performance.md) | 0.5708 | 1.03 | 1.18 | 4.54 | 2.32 | 3.21 | 10.29 | 15.69 | 17.11 |
| [asin](functions/asin/Performance.md) | 0.8320 | 0.87 | 0.88 | 4.38 | 1.73 | 1.70 | 13.85 | — | 9.45 |
| [acos](functions/acos/Performance.md) | 0.8174 | 0.97 | 0.94 | 4.56 | 1.86 | 1.94 | 14.44 | — | 9.84 |
| [atan](functions/atan/Performance.md) | 0.5640 | 1.00 | 1.07 | 5.77 | 3.35 | 2.21 | 16.48 | — | 8.80 |
| [atan2](functions/atan2/Performance.md) | 0.7769 | 1.29 | 1.92 | 5.36 | 4.03 | 2.70 | 13.78 | — | 19.22 |
| [sinh](functions/sinh/Performance.md) | 0.8696 | 0.79 | 0.85 | 1.69 | 2.15 | 2.46 | 5.08 | — | 8.49 |
| [cosh](functions/cosh/Performance.md) | 0.7510 | 0.56 | 0.62 | 1.13 | 1.59 | 2.44 | 3.60 | — | 7.74 |
| [tanh](functions/tanh/Performance.md) | 0.6436 | 0.97 | 1.63 | 3.10 | 2.80 | 3.18 | 6.97 | — | 17.02 |
| [asinh](functions/asinh/Performance.md) | 0.7383 | 1.52 | 1.50 | 8.06 | 5.12 | 5.18 | 29.55 | 25.52 | 11.09 |
| [acosh](functions/acosh/Performance.md) | 0.7349 | 1.65 | 1.69 | 8.60 | 4.59 | 5.51 | 30.44 | 23.70 | 10.02 |
| [atanh](functions/atanh/Performance.md) | 0.7178 | 1.29 | 1.20 | 5.42 | 4.30 | 2.92 | 19.15 | 15.95 | 19.21 |

¹ Highway hp algorithm (max 0.505 ulp on AVX-512 with an 8-entry table, 0.511 ulp elsewhere), not Intel's; Intel's
HA exp has max 0.80 ulp. The C translation of Intel's algorithm runs at 0.43 ns on AVX-512.
² Intel's shipped libsvml AVX-512 pow is a different build from the open-source code ported here (it gives
different results); the AVX2 column compares like with like more closely.

## Known gaps

- AVX-512: the emulations of `vreducepd`, `vscalefpd` and rounding-mode arithmetic still run on AVX-512 targets
  (exp2, expm1, pow, cbrt, the log family), where the instructions themselves would be cheaper.
- pow is latency-bound (long dependency chains of exact double-double steps): 2.8× Intel's AVX-512 time and
  2.5× its AVX2 time.
- AVX2: the exact emulation of `vrcp14pd`/`vrsqrt14pd` makes asin, acos, atan, atan2 and the inverse
  hyperbolic functions 1.3–2.6× slower than Intel's (differently designed) AVX2 code.
- Large trig arguments (|x| > 2^16): the Highway large-argument path is about 2× slower than Intel's.

Local machine: Intel Core i5-9500 (Coffee Lake, AVX2), see the per-function pages.
