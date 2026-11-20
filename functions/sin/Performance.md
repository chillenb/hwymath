# sin (double precision): performance

Accuracy and local throughput of the sin implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (real AVX-512, AArch64) are
added at integration; timings under Intel SDE are meaningless and not reported.

## How it is measured

- **Accuracy**: `make sin-hwy` (`tests/hwy_test.cc`), errors in ulps of the correctly rounded result against
  `sinl` (x86-64 `long double`, 64-bit mantissa, which slightly undercounts results that are not
  correctly rounded), over the 20,035,562 inputs of `functions/sin/tests/inputs.h`, split into the
  main path (|x| <= 2^24, plus +-inf/NaN, which must give exactly the reference's NaN) and the
  large-argument path (finite |x| > 2^24). "not CR": error > 0.5 ulp.
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
| Highway AVX2 | SVML (bit-identical) | 4 | 0.5200 | 0.0619 % | 0 | 0.5024 | 0.0029 % | 0 | 2.17 | 17.9 |
| Highway SSE4 / SSSE3 / SSE2 | `std::sin` | 2 | 0.5151 | 0.0727 % | 0 | 0.5151 | 0.1189 % | 0 | 15.08 | 63.0 |
| Intel SVML HA AVX2 |  |  | 0.5029 | 0.0023 % | 0 | 0.5024 | 0.0029 % | 0 | 3.10 | 9.1 |
| SLEEF u10 AVX2 |  |  | 0.7539 | 1.0226 % | 0 | 0.8564 | 1.9078 % | 0 | 4.46 | 18.8 |
| libmvec AVX2 |  |  | 3.2896 | 31.8505 % | 1785937 | 0.5151 | 0.1189 % | 0 | 1.07 | 70.3 |
| glibc sin (scalar) |  | 1 | 0.5151 | 0.0727 % | 0 | 0.5151 | 0.1189 % | 0 | 14.99 | 63.2 |

Worst inputs of the SVML algorithm: x = `0x1.7ca4bd587376p+21` (main path), `0x1.83c15911b9b7p+26` (large path).
All targets running the SVML algorithm (AVX2 here; AVX3, AVX3_DL, AVX3_SPR under SDE) give identical
results, which are those of the C translation and of Intel's assembly, so their accuracy is the same;
Intel's libsvml AVX-512 version (`__svml_sin8_ha_z0`, under SDE) agrees with them on all inputs except
sin(+-inf) (NaN sign).

The large-argument path of the Highway version (~18 ns/value) loads its tables (3 words of 2/pi per
exponent, 4 words per 512-entry row) with scalar loads instead of gathers; Intel's AVX2 libsvml takes
~8.5 ns there.
