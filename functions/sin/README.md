# sin (double precision)

Source studied: [`svml_z0_sin_d_ha.s`](https://github.com/numpy/SVML/blob/9a43e744b15006d63e8d1b57626dfbe6b984f5ee/linux/avx512/svml_z0_sin_d_ha.s)
(AVX-512, high accuracy). Accuracy and speed: [Performance.md](Performance.md). cos
([`functions/cos`](../cos/README.md)) has the same structure and tan ([`functions/tan`](../tan/README.md))
the same large-argument reduction; the code they share lives here.

| path | contents |
|---|---|
| `svml/svml_sin_d_ha.c` | bit-identical C/AVX-512 translation of the assembly and its scalar callout |
| `svml/svml_sin_d_ha_constants.inc` | its constants, generated from the assembly by `tools/gen_svml_constants.py sin` |
| `svml/svml_trig_reduction.h` | **shared by sin, cos, tan**: the Payne-Hanek large-argument reduction (C/AVX-512) and the evaluation tail of the sin/cos large path |
| `svml/svml_trig_reduction_data.inc` | **shared**: its 2048x192-bit table of 2/pi and constants (`gen_svml_constants.py reduction`; identical in all three sources) |
| `hwy/hwy_sin.cc` | Highway version |
| `hwy/trig_hwy-inl.h`, `hwy/trig_tables.h` | **shared**: Highway table lookups (vpermt2pd emulation) and the large-argument path |
| `tests/` | `test_rare.c`, `test_sin8.c`, `hwy_test.cc`, `extras.c`; **shared** templates `*_tmpl.h` and `inputs.h` |
| `tools/gen_svml_constants.py` | constant generator for sin, cos and the reduction data |
| `trig_common.mk` | **shared** make variables (included by the sin, cos and tan `function.mk`) |
| `sin_d.h` | declarations of the C translation |

**Why three directories with shared files in `sin/`:** sin and cos share most of their code (the same
tables, the same large path; only the main-path formula and its operation order differ), tan shares the
large-argument reduction only. Each function keeps its own translation (the operation order differs
instruction by instruction), its own Highway kernel, tests and documentation; the reduction, the Highway
building blocks, the test templates and the constant generator exist once, here.

## Algorithm

From the assembly's header comment and the code:

**Main path** (|x| <= 2^24):
- `y = x*16/pi + 1.5*2^52` (FMA, round to nearest): N = round(x*16/pi) sits in the low bits of y, the
  table index j = N mod 16 (vpermt2pd uses the low 4 bits of y directly).
- R + Rlo = x - N*(PI1 + PI2 + PI3) with pi/16 split in three parts (|R| <= pi/32); Rlo collects the
  rounding errors of the two-step reduction (exact FMA error terms).
- With B = N*pi/16: Sh + Sl = sin B, Ch + Cl = cos B from three 16-entry tables (the low parts share one
  64-bit word: Cl is the word read as a double, Sl its low half shifted up by 32 bits).
- sin(B + R) = Sh + R*Ch [with its exact rounding error] + R^2*(Sh*Cp + R*Ch*Sp) + Sl
  &plus; Rlo*(Ch - Sh*R) + R*Cl, where cos R - 1 = R^2*Cp(R^2) and sin R - R = R^3*Sp(R^2) (degree-3 polynomials in R^2).
- Bit 4 of N (B in [pi, 2pi) mod 2pi) flips the sign. For x = -0 the code flips bit 4 of N once more
  (integer compare with -0.0), so that sin(-0) = -0.

**Large path** (|x| > 2^24, +-inf, NaN; if any lane needs it, all 8 lanes are recomputed and merged):
a Payne-Hanek reduction modulo pi/256 (`svml/svml_trig_reduction.h`):
- e = biased exponent of |x|, M = significand as a 53-bit integer. Three 64-bit words T0..T2 of a table
  indexed by e (2048 entries, 48 KiB, read with `vgatherdpd`) hold the 192 bits of 1/(2pi) that matter
  for this exponent (the bits that would only contribute integer multiples of 2pi are left out).
- The product M*T (truncated: ten 32x32->64-bit `vpmullq` of the 32-bit halves, carries propagated by
  hand) gives the fraction of x/(2pi) as a 128-bit fixed-point number hi:lo.
- Bits 127..76 become the double F = 1.f (implicit 1 added through the exponent bits); adding
  1.5*2^43 rounds F to a multiple of 2^-9: N = that multiple mod 512 (x = N*2pi/512 + r), and the
  remainder F - N/512 is exact. The next 52 and 24 bits become two more doubles (OR the exponent of
  2^-52 / 2^-104 into the bits, subtract the same power of 2). Summing the three pieces into RHi + RLo and
  multiplying by 2pi (two parts, FMA for the exact product error) gives r + rlo, |r| <= pi/512.
- Inputs with |x| <= 2^-20 would pass through unreduced (dead code here, since |x| > 2^24).
- Evaluation with a 512-entry table of {T0, T1, T2, T3} per N: sin B = T2 + T3, cos B = T0 + T1 with T1 a
  power of 2 (2^round(log2 cos B)), and f(B + r) = T2 + (T0 + T1)*r + T2*(cos r - 1) + (T0 + T1)*(sin r - r)
  &plus; T3 + rlo*(T0 + T1 - T2*r), with an exact split of T2 + T1*r + T0*r and degree-2/3 polynomials. The
  table is read with four gathers; the sign of x is applied at the end.

**Callout** (`__svml_dsin_ha_cout_rare_internal`): called for the lanes with x = +-inf or NaN (the mask
is computed in the large path). sin(+-inf) = 0*x (NaN, invalid flag, return code 1); sin(NaN) = x*x;
other inputs: nothing written. The vector result of these lanes is replaced.

## The C translation (`svml/svml_sin_d_ha.c`)

`svml_sin8_ha(__m512d)` translates the vector routine with AVX-512F/DQ intrinsics: same algorithm,
constants (hex-float literals generated from the assembly), operation order and rounding mode of every
instruction (`_mm512_*_round_pd` with `{rn-sae}`), so its results are bit-identical to Intel's
(`test_sin8`). `svml_sin_ha_array(x, y, n)`: loop with a masked tail. `svml_dsin_ha_cout_rare()`: the
callout in portable C.

Differences from the assembly that do not affect results:
- bitwise operations use `vpandq`/`vpxorq` where the assembly sometimes uses `vandpd`/`vxorpd`;
- the large path is a separate `noinline` function (in the assembly it is a block after the `ret`);
- the order of operands of commutative operations is the compiler's: it could only matter for NaN
  payloads, and NaN results only occur in lanes the callout replaces;
- the callout stores the x86 default NaN `0xfff8000000000000` for sin(+-inf) explicitly (after computing
  0*inf at run time for the invalid flag), so that it gives the same bits on AArch64, whose default NaN
  is positive.

Intel's libsvml (`__svml_sin8_ha_z0`, conda-forge `intel-cmplr-lib-rt`) agrees with the translation on
all test inputs except sin(+-inf), where it returns +NaN instead of the assembly's -NaN.

## Highway version (`hwy/hwy_sin.cc`)

The same algorithm, operation for operation, with [Google Highway](https://github.com/google/highway)
1.4.0 and run-time dispatch:

- targets with true FMA (`HWY_NATIVE_FMA`: AVX2, AVX3*, NEON, SVE*): bit-identical to Intel's AVX-512
  code;
- targets without FMA (SSE2/SSSE3/SSE4, EMU128): `std::sin` (the FMA-based error-free steps cannot be
  reproduced).

Emulation of each AVX-512 feature (all exact):

| AVX-512 feature | emulation |
|---|---|
| `{rn-sae}` arithmetic (add, sub, mul, FMA variants) | Highway's `Add`, `Sub`, `Mul`, `MulAdd`, `NegMulAdd`, `MulSub`, `NegMulSub` (round to nearest; compiled with `-ffp-contract=off`) |
| `vpermt2pd` (16-entry tables, index = low 4 bits) | `Lookup<16>` (`hwy/trig_hwy-inl.h`): `TableLookupLanes` if one vector holds the table, else a tree of `TwoTablesLookupLanes` and sign-bit selects (AVX-512: one `vpermt2pd`; AVX2: 2 x 2 `vpermd` + blends; 128-bit: 4 `TBL2`/shuffles + 3 selects); scalable vectors shorter than 16 lanes: `GatherIndex` (SVE) |
| masked `vxorpd` (sin(-0)), `vpandnq` mask tricks, merges with `k1` | `IfThenElse`, `IfThenElseZero`, `IfThenZeroElse` |
| `vcmppd` with unordered predicates (`NLE_UQ`, `NGT_UQ`) | `Not(Le(..))`, `Not(Gt(..))` (true for NaN, like the predicates) |
| `vpmullq` of 32-bit values | `MulEven` on the 32-bit halves (exact 64-bit products) |
| `vgatherdpd` (reduction table by exponent, 512-entry table) | scalar loads into an aligned buffer (`LoadIndexed`): only in the rare large path; no gather instructions (slow with the Gather Data Sampling mitigation) |
| integer shifts/adds/ands of the reduction | the same Highway integer operations on 64-bit lanes |

Vectors with a lane |x| > 2^24 (or inf/NaN) are recomputed by an out-of-line function with the large
path; lanes with inf/NaN then call the C callout. **Deviations: none** on FMA targets. Bit-identity was
verified under Intel SDE (`-spr`): AVX3_SPR, AVX3_DL, AVX3 and AVX2 against the C translation on
2,035,562 inputs (`make sin-avx512-test`), and natively the AVX2 hash equals the hash of the
C translation (`make sin-hashes`). AArch64 NEON/SVE are expected to agree (same IEEE operations; the
callout gives the same NaN bits on both architectures) but were not run here.

## Tests (`make sin-test`, `make sin-avx512-test RUN="$SDE -spr --"`)

| test | compares with | result |
|---|---|---|
| `test_rare` | Intel's scalar `__svml_dsin_ha_cout_rare_internal`: results, whether written, return codes, FP exception flags | 20,035,562 inputs, 0 mismatches |
| `test_sin8` | Intel's `__svml_sin8_ha`, all 8 lanes bit for bit (needs AVX-512F/DQ: Intel SDE here) | 20,035,568 inputs (7,350,780 in the large path), 0 mismatches |
| `hwy_test` | every Highway target: error < 1 ulp, exact special values, identical results on all SVML-algorithm targets, bit-identity with the C translation when AVX-512 is available; comparisons with Intel libsvml, SLEEF, libmvec, glibc | all checks pass |
| `make sin-insns` | arithmetic instruction mix of the compiled main path vs the assembly | identical |

Inputs (`tests/inputs.h`, shared with cos and tan): special values (+-0, +-inf, quiet and signaling NaNs,
denormals, DBL_MAX), +-50 doubles around every threshold (2^24, tan's 65536.0035, 2^-20, 2^52/53/63/64,
2^1023, ...), the doubles next to k*pi/2 (+-2 ulps) for |k| <= 2000 and for 2000 random k up to 2^49,
the main-path rounding boundaries (n + 1/2)*pi/16, the large-path rounding boundaries (2M + 1)*pi/512
near 2^17 and 2^25, 6381956970095103*2^797 (the double closest to a multiple of pi/2), and random values:
uniform in [-10, 10], [-1e5, 1e5], [-2^24, 2^24], [-2^30, 2^30], log-uniform over all exponents, and
arbitrary bit patterns.
