# exp2 (double precision): performance

Accuracy and throughput of the exp2 implementations in this directory and of the libraries they are
compared with. Descriptions of the versions: [README.md](README.md). CI numbers (AVX-512, AArch64)
are added at integration. Timings under Intel SDE are meaningless and are not reported.

## How it is measured

- **Accuracy**: errors in ulps of the correctly rounded result, against `exp2l` (`long double`), over
  the test inputs of `tests/inputs.h` (20,016,321 inputs: special values, every threshold with
  neighbouring doubles, table boundaries, uniform in [−1100, 1100] and [−1, 1], random bit patterns).
  The statistics cover the inputs with a finite, nonzero result. "not CR" means not correctly rounded
  (error > 0.5 ulp). "special" counts NaN/inf/zero mismatches. The x86-64 `long double` has a 64-bit
  mantissa, so the x86-64 "not CR" percentages are slightly optimistic. With AArch64's binary128
  (QEMU, `make exp2-xcheck`, 2M inputs), the SVML algorithm has max 0.7477 ulp and 3.5569 % not CR.
- **Throughput**: ns/value for 4096 random doubles in [−1000, 1000] (cache-resident, no callouts),
  best of 7 runs of 2000 repetitions, array functions (`make exp2-hwy`).
- **Toolchain**: conda-forge GCC via pixi, Highway 1.4, glibc 2.28 sysroot.

## Core i5-9500 (AVX2), local

Accuracy: `hwy_test --no-bench` (20M random inputs). Timings: `hwy_test 2000000` on a shared machine,
so expect some noise.

| implementation | algorithm | lanes | max ulp | not CR | > 1 ulp | special | ns/value |
|---|---|---|---|---|---|---|---|
| **Highway AVX2** | SVML (bit-identical) | 4 | **0.7520** | 3.5523 % | 0 | 0 | **1.72** |
| Highway SSE4/SSSE3/SSE2 | `std::exp2` | 2 | 0.5063 | 0.0306 % | 0 | 0 | 8.30 |
| Highway contrib `Exp2`, AVX2 (x ≤ 1024) | | | 1.0234 | 6.2059 % | 13 | 4905 | 1.23 |
| Intel SVML HA AVX2 (`__svml_exp24_ha_l9`) | | | 0.7520 | 0.3317 % | 0 | 0 | 1.01 |
| SLEEF u10 AVX2 | | | 0.8481 | 4.4560 % | 0 | 0 | 1.30 |
| glibc `exp2` (scalar) | | 1 | 0.5063 | 0.0306 % | 0 | 0 | 8.30 |

## AVX-512 (accuracy only, Intel SDE `-spr`, 2,016,321 inputs)

All Highway FMA targets (AVX3_SPR, AVX3_DL, AVX3, AVX2) give results bit-identical to the C/AVX-512
translation, and so to `__svml_exp28_ha`: max 0.7476 ulp, 3.5366 % not CR, none > 1 ulp. On the same
inputs, Intel's libsvml `__svml_exp28_ha_z0` (a different build) has max 0.7629 ulp and 3.5400 % not CR.
SLEEF u10 AVX-512 has max 0.8154 ulp.

## Error by path (Highway AVX2, 20M inputs)

| lanes | max ulp | not CR |
|---|---|---|
| vector path, \|x\| < 1 | 0.7251 | 2.79 % |
| vector path, 1 ≤ \|x\| < 1022 | 0.6943 | 5.25 % |
| callout, normal results | 0.5137 | 0.59 % |
| callout, denormal results (x < −1022; double rounding) | 0.7520 | 0.94 % |

The vector path's 3.5 % not-CR rate comes from the algorithm. The small part Tl + (Th·r)·P(r), up to
about 0.045·Th, is computed in plain double arithmetic (Th·r and the Horner steps each round), and it
is added to Th with one final rounding. Intel's AVX2 code (`_l9`) is a different algorithm and is
correctly rounded far more often.
