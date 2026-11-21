# pow: accuracy and speed

Measured on the development machine: Intel Core i5-9500 (Coffee Lake, AVX2/FMA, no AVX-512), shared with
other jobs (load ~3–4 on 6 cores), so timings are only indicative. AVX-512 rows come from Intel SDE 10.13
(`-spr`); SDE timings are meaningless and are not shown. CI numbers on real hardware are added at
integration.

## Bit-identity

| check | inputs | result |
|---|---|---|
| `test_rare`: C callout vs `__svml_dpow_ha_cout_rare_internal` (results, return codes, FP flags) | 20,095,065 pairs | 0 mismatches |
| `test_pow8`: `svml_pow8_ha` vs `__svml_pow8_ha` (SDE) | 2,095,072 pairs | 0 mismatches |
| `hwy_test` under SDE: Highway AVX3_SPR, AVX3_DL, AVX3, AVX2 vs the C translation | 2,095,065 pairs | 0 differences on every target |
| emulated AVX-512 operations vs the instructions (SDE), on every SVML target | 0.66 M rcp, 0.4 M each for mul/fma/add RZ, reduce and add RD | 0 mismatches |
| `make pow-hashes` natively (AVX2) vs C translation under SDE | 2,095,065 pairs | same hash `95289bd0174d4859` |
| `make pow-xcheck`: AArch64 NEON, SVE and SVE2 (128/256/512/2048-bit) under QEMU vs x86-64 AVX2 | 395,065 pairs | bit-identical (NaNs canonical, see below) |

The only difference between architectures is the sign of the default NaN that 0/0 produces for x < 0 with
non-integer y: negative on x86-64 (as in Intel's code), positive on AArch64. Hashes and xcheck outputs
therefore map all NaNs to one canonical NaN.

## Accuracy

The reference is `powl` (x87 80-bit). 2,095,065 input pairs from `tests/inputs.h`, excluding signaling-NaN
inputs. Results that `powl` rounds to 0 or ±inf must match exactly ("special" column). Errors near
overflow are measured against `powl`'s own result, so its 64-bit mantissa slightly undercounts results
that are not correctly rounded.

| implementation | max ulp | not correctly rounded | > 1 ulp | special mismatches | worst case |
|---|---|---|---|---|---|
| Highway, all SVML targets (AVX2, AVX3, AVX3_DL, AVX3_SPR) = C translation = Intel's assembly | 0.7490 | 2.7918 % | 0 | 0 | pow(1 + 2^−52, −1.385·2^61): denormal result, scalar callout |
| Highway SSE2/SSSE3/SSE4 (`std::pow` fallback) | 0.5068 | 0.0359 % | 0 | 0 | |
| Intel libsvml HA AVX2 (`__svml_pow4_ha_l9`) | 0.7490 | 0.1428 % | 0 | 0 | same as above |
| Intel libsvml HA AVX-512 (`__svml_pow8_ha_z0`, SDE; not the same build as the studied assembly) | 0.8247 | 2.8303 % | 0 | 0 | pow(2, −1022.0…) |
| SLEEF u10 (AVX2 and AVX-512) | 0.7749 | 1.0917 % | 0 | 0 | |
| glibc libmvec (AVX2 and AVX-512) | 0.9961 | 10.7164 % | 0 | 0 | |
| glibc `pow` (scalar) | 0.5068 | 0.0359 % | 0 | 0 | |

On the vector path, the error stays below 0.75 ulp. The worst cases are denormal results, which go through
the scalar callout.

## Speed (AVX2, i5-9500, ns per value, best of 7)

4096 pairs with x in [0.01, 100] and y in [−50, 50], so no lane needs the callout:

| implementation | ns/value |
|---|---|
| Highway AVX2 (bit-identical to Intel's AVX-512 code) | 14.3 |
| Intel libsvml HA AVX2 (`__svml_pow4_ha_l9`, different algorithm) | 7.0 |
| SLEEF u10 AVX2 | 11.0 |
| glibc libmvec AVX2 | 10.2 |
| glibc `pow`, scalar (also the SSE fallback) | 13.5 |

The AVX2 version pays for exactness. Each 32-entry table lookup costs 8 `vpermd` and 7 blends, and the
four round-toward-zero emulations add about 40 operations, mostly the ErrFma error-sign computation for
the RZ fma. The dependency chain is about 230 cycles long. `perf` shows the reservation station full for
half of the cycles while the execution ports are only ~37 % busy: the kernel is latency-bound.
Unrolling by two vectors (with GCC's pre-RA scheduling) did not help. On AVX-512 the lookups are the
same single `vpermt2pd` + blend as in the assembly, and the other AVX-512-only operations are the
instructions themselves (README: the AVX3 loop has the C translation's 67 vector arithmetic instructions).
