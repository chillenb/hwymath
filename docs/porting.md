# Porting an SVML high-accuracy function

How each function in `functions/` is studied. `functions/exp/` is the reference implementation of this
guide: read its README, `function.mk`, `svml/svml_exp_d_ha.c`, `hwy/hwy_exp.cc` and `tests/` first.

## Deliverables (double precision)

For a function `f` with source `third_party/SVML/linux/avx512/svml_z0_f_d_ha.s` (vector entry
`__svml_f8_ha`, scalar callout `__svml_df_ha_cout_rare_internal`, data `__svml_df_ha_data_internal*`):

1. **Bit-identical C translation** `functions/f/svml/svml_f_d_ha.c`: AVX-512 intrinsics, same algorithm,
   constants, operation order and per-instruction rounding/SAE (`_mm512_*_round_*`), plus the scalar
   callout translated to portable C (`svml_df_ha_cout_rare`, compiled without FMA contraction; values
   that must raise FP flags at run time need `volatile`). Only the vector function gets
   `__attribute__((target("avx512f[,avx512dq]")))`, so the file builds and the callout runs anywhere;
   guard x86-only parts with `#if defined(__x86_64__)`. Constants as hex-float literals generated from
   the assembly by `functions/f/tools/gen_svml_constants.py` (use `common/tools/extract_tables.py`)
   into `functions/f/svml/svml_f_d_ha_constants.inc` (committed). Header `functions/f/f_d.h`.
2. **Highway version** `functions/f/hwy/hwy_f.cc` (dynamic dispatch; see `functions/exp/hwy/`): a port of
   the same SVML algorithm. Goal: bit-identical to Intel's AVX-512 code on every target with true FMA
   (`HWY_NATIVE_FMA`), by emulating AVX-512-only operations exactly — rounding modes via an FMA plus an
   exact error-sign correction, `vscalefpd` via exponent bits, `vgetexp/vgetmant` via integer ops (care
   with denormals), `vpermt2pd` via `TableLookupLanes`/`TwoTablesLookupLanes` (tables must be
   `constexpr`), `vfixupimm`/`vfpclass`/`vrange` via compares and selects, `vrcp14pd`/`vrsqrt14pd` with the
   shared exact emulation (`common/include/hwymath/hwy/avx512_emul-inl.h`, after Intel's reference
   `recip14.c`; see `common/README.md`). On AVX-512 targets (`#if HWY_TARGET <= HWY_AVX3`) use the
   instruction itself wherever Highway or an intrinsic on `.raw` provides it: it is exact by definition
   and much cheaper than the emulation. Should exactness really be impossible somewhere, keep the
   algorithm, document the deviation, and require the documented accuracy (max error measured, < 1 ulp
   for HA). Targets without FMA fall back to the `std::`
   function. No gathers (very slow with the GDS microcode mitigation). Lanes needing the scalar callout use
   the C callout from step 1. Export `extern "C"` array functions; call `HWY_DYNAMIC_DISPATCH` inside the
   namespace with unqualified names. Two-argument functions (atan2, pow) take two input arrays.
3. **Tests** in `functions/f/tests/` (exit status != 0 on failure; 77 = skipped):
   - `inputs.h`: special values (±0, ±inf, NaNs, domain edges, every threshold of the algorithm and of the
     callout with neighbouring doubles), table boundaries, random values over the useful range and random
     bit patterns. Deterministic (`common/include/hwymath/test_util.h`).
   - `test_rare.c`: C callout vs the assembly's callout (made global with `sed` as in exp): results,
     return codes and FP exception flags, bit for bit.
   - `test_f8.c`: AVX-512 translation vs `__svml_f8_ha` bit for bit (skip with 77 without AVX-512F).
   - `hwy_test.cc`: every Highway target that is compiled in and supported: accuracy vs the `long double`
     reference (`common/include/hwymath/accuracy.h`), bit-identity between targets, bit-identity with the C
     translation where claimed (AVX-512 targets under SDE; AVX2/NEON/SVE via hashes), special values,
     `--no-bench`, and `--hash FILE` writing `arch f target variant n hash` lines (and
     `arch f inputs - n hash`) for the CI cross-machine check. Optional comparisons: Intel libsvml
     (`__svml_f4_ha_l9`, `__svml_f8_ha_z0`; only call it, never disassemble it), SLEEF, glibc libmvec.
4. **`functions/f/function.mk`**: modelled on exp's: `build/f/`, targets `f-all f-test f-hwy f-hashes`,
   and on x86-64 `f-insns` (`common/tools/compare_insns.sh`) and `f-avx512-test` (honouring `$(RUN)`),
   appended to the aggregate variables. The top-level Makefile includes it automatically.
5. **`functions/f/README.md`** (algorithm as understood from the assembly, versions, emulation of each
   AVX-512 feature, deviations, tests) and **`functions/f/Performance.md`** (accuracy measured on this
   machine; local AVX2 timings; SDE timings are meaningless — CI numbers are added at integration).

Several closely related functions may share one directory only if they share most code; otherwise one
directory per function.

## Shared code

Use what is in `common/` (overview: `common/README.md`) instead of copying it: `accuracy.h` (ulp errors,
statistics, hashes, timing), `hwy/avx512_emul-inl.h` (`vrcp14pd`, `vrsqrt14pd`, `vgetexppd`, 16-entry
`vpermt2pd`), `rcp14_ref.h` (scalar `vrcp14pd` reference), `tools/xcheck.sh` (QEMU cross-check; the
function's `tests/hwy_xcheck.cc` writes the inputs and results). Code shared only within a family lives in
the family's first directory (`functions/log/`, `functions/sin/`).

## Verification before calling a function done

```bash
make f-test                                   # translation callout, Highway targets on this CPU
make f-avx512-test RUN="$SDE -spr --"         # AVX-512 translation vs Intel's assembly, Highway AVX-512 targets
make f-insns                                  # instruction mix of the translation vs the assembly
make f-hashes HASH_FILE=/tmp/...              # AVX2 hash must equal the translation's results if bit-identical
```

`SDE=/home/claude/opt/sde-external-10.13.1-2026-07-28-lin/sde64` (Intel SDE; `-spr` = Sapphire Rapids).
Keep SDE runs to ~2M inputs. Run all of the above clean from `make clean`-state before reporting.

## Environment and rules for parallel work

- Each function family is worked on in its own git worktree and branch. Activate the shared pixi
  environment instead of installing a new one:
  `eval "$(~/.pixi/bin/pixi shell-hook --manifest-path /home/claude/hwymath/pixi.toml)"`
  (then plain `make ...` in the worktree).
- The machine (6 cores, 15 GB) is shared by several agents: use at most `make -j2`, prefix long test
  runs with `nice`, and keep random input counts moderate.
- Only create or modify files under `functions/<your functions>/`. Do not modify `common/`, the
  top-level `Makefile`, `README.md`, `Performance.md`, `docs/`, `.github/` or `pixi.*`; if you need a
  shared change, put a local version under your function directory and describe the proposed shared
  change in your report.
- Commit on your branch (author identity is configured); do not push.

## Gotchas already hit

- Compile with `-ffp-contract=off` (set in `common/mk/config.mk`).
- The SVML constant macros (`vC1`, ...) collide with identifiers inside Highway's headers: turn the needed
  ones into `constexpr` in a header and `#undef` all macros before including Highway.
- Highway's static AVX2 target also needs AES/PCLMUL (`-march=native` for static x86 builds).
- `sleef.h` only declares AVX2/AVX-512 functions when compiled with those ISAs: declare them manually.
- The x86-64 `long double` reference (64-bit mantissa) slightly undercounts not-correctly-rounded results.
- The conda sysroot is glibc 2.28: scalar libm functions may bind to older, slower symbol versions.
- Test inputs must not depend on the architecture: never call the random generator twice in one expression
  (`f(rng(), rng())`: the evaluation order is unspecified and differs between GCC on x86-64 and AArch64),
  and keep `-ffp-contract=off` for the input generator.
- The default NaN of an invalid operation is negative on x86-64 and positive on AArch64; hash results with
  `acc_hash_nan_canonical` when the function creates NaNs from non-NaN inputs (pow), and ignore NaN payloads
  when comparing with Intel's code.
- Highway's `NegMulSub(a, b, c)` on NEON is `Neg(MulAdd(a, b, c))`, which gives the opposite signed zero to
  x86 and SVE when the result is zero (tan(−0) came out +0 on NEON). Use `NegMulAdd(a, b, Neg(c))` where the
  sign of a zero result matters.
- SVE vectors are sizeless: no arrays of vectors and no vectors in structs (use template recursion or
  separate variables); compile the AArch64 multi-target build (`make <f>-xcheck`, or CI) before pushing.
- `compare_insns.sh` compares instruction mixes, so the C translation must compile to the same instructions:
  empty `asm` barriers stop GCC from rewriting sign manipulations (`x ^ |x|`) as `and`/`vpternlogq` or
  replacing loads of all-zero table halves.
- On AVX-512, the emulations of `vreducepd`, `vscalefpd` and rounding-mode arithmetic (`{rz-sae}`,
  `{rd-sae}`) cost several instructions each, and AVX-512 code is usually limited by the two vector ports:
  instruction count is time. Prefer the instructions on AVX-512 targets.
