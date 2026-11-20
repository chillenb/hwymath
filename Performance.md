# Performance

Accuracy and speed are documented per function, with the measurement method and the machines:

| function | details | best version here | max error | Intel SVML HA (widest) | SLEEF u10 (widest) |
|---|---|---|---|---|---|
| exp (double) | [functions/exp/Performance.md](functions/exp/Performance.md) | Highway hp | 0.505 ulp (AVX-512), 0.511 ulp (AVX2, NEON, SVE) | 0.800 ulp (AVX-512) | 0.877 ulp |

Machines used so far: Intel Xeon 6975P-C (Granite Rapids, AVX-512) and Neoverse-N2 (SVE2, 128-bit) in
GitHub Actions; Intel Core i5-9500 (Coffee Lake, AVX2) locally.
