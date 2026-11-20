# tan (double precision): performance

Accuracy and local throughput of the tan implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (real AVX-512, AArch64) are
added at integration; timings under Intel SDE are meaningless and not reported.

## How it is measured

- **Accuracy**: `make tan-hwy` (`tests/hwy_test.cc`), errors in ulps of the correctly rounded result against
  `tanl` (x86-64 `long double`, 64-bit mantissa, which slightly undercounts results that are not
  correctly rounded), over the 20,035,562 inputs of `functions/sin/tests/inputs.h`, split into the
  main path (|x| <= 65536.0035, plus +-inf/NaN, which must give exactly the reference's NaN) and the
  large-argument path (finite |x| > 65536.0035). "not CR": error > 0.5 ulp.
- **Throughput**: ns/value for 4096 cache-resident inputs, uniform in [-100, 100] ("ns") and log-uniform
  in [2^30, 2^1000] ("ns huge"), best of 7 runs. Array functions; Intel's libsvml, SLEEF and libmvec
  are called once per vector.
- **Machine**: Intel Core i5-9500 (Coffee Lake, AVX2/FMA, no AVX-512), shared with other jobs (expect
  some noise), conda-forge GCC 15 via pixi with the glibc 2.28 sysroot: the scalar glibc functions and
  libmvec are old versions (the large errors of scalar cos/tan for huge arguments below come from the
  older symbol version; the system's current glibc rounds the worst of them, 6381956970095103*2^797,
  correctly).

## Results (Core i5-9500)

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | max ulp, huge | not CR, huge | > 1 ulp, huge | ns | ns huge |
|---|---|---|---|---|---|---|---|---|---|---|
| Highway AVX2 | SVML (bit-identical) | 4 | 0.5708 | 0.1753 % | 0 | 0.5703 | 0.3759 % | 0 | 5.06 | 23.4 |
| Highway SSE4 / SSSE3 / SSE2 | `std::tan` | 2 | 0.6821 | 0.1054 % | 0 | 14.3608 | 0.2277 % | 3 | 19.74 | 70.1 |
| Intel SVML HA AVX2 |  |  | 0.5156 | 0.0727 % | 0 | 0.5210 | 0.1564 % | 0 | 3.29 | 8.8 |
| SLEEF u10 AVX2 |  |  | 1.0000 | 0.5705 % | 0 | 0.6255 | 1.1133 % | 0 | 7.29 | 18.1 |
| glibc tan (scalar) |  | 1 | 0.6821 | 0.1054 % | 0 | 14.3608 | 0.2277 % | 3 | 19.24 | 68.3 |

Worst inputs of the SVML algorithm: x = `0x1.3c546555ea0d4p+3` (main path), `0x1.c1471f9502fe8p+22` (large path).
All targets running the SVML algorithm (AVX2 here; AVX3, AVX3_DL, AVX3_SPR under SDE) give identical
results, which are those of the C translation and of Intel's assembly, so their accuracy is the same;
Intel's libsvml AVX-512 version (`__svml_tan8_ha_z0`, under SDE) agrees with them on all inputs except
tan(+-inf) (NaN sign).

SLEEF's "1.0000" is an error just below 1 ulp (no result above 1 ulp).

The Highway AVX2 version is slower than Intel's AVX2 libsvml (a different, AVX2-specific algorithm with
a division) mainly because of the exact `vrcp14pd` emulation: replacing it by a division brings the
same code to ~2.6 ns/value (but changes 0.2 % of the results; see README.md). The large-argument path
(~23 ns/value) loads its tables with scalar loads instead of gathers.
