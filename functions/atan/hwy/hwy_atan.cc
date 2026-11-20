// Portable (Google Highway) version of SVML's high-accuracy atan (svml/svml_atan_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation by operation, with exact
//     equivalents of the AVX-512 features: vreducepd (a - round(4a)/4) = a - (RN(a + Shifter) - Shifter)
//     (the same rounding to a multiple of 1/4, ties to even; the subtraction is exact), the two-source
//     permutes of the 32-entry tables = TableLookupLanes/TwoTablesLookupLanes (vector length dependent:
//     one vector, two vectors, or several two-vector lookups merged with selects, e.g. four on AVX2),
//     vminpd = Min, masked operations = IfThenElse. vrcp14pd, whose result bits are implementation-
//     defined, is replaced by trunc16(RN(1/Y)) (16 significant bits, like Intel's, relative error
//     < 2^-16 + 2^-53 instead of 2^-14), so the results are identical on all these targets and machines
//     (= svml_atan8_ha_portable), and differ from Intel's in the last bit for a small fraction of the
//     inputs (see README.md);
//     AtanNativeArray uses vrcp14pd itself on AVX-512 targets (bit-identical to Intel's code on the same
//     CPU) and is the same as AtanArray elsewhere;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::atan.
// No scalar callout; NaN lanes return the quieted input, like the SVML code.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "functions/atan/hwy/hwy_atan_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/atan/hwy/hwy_atan.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwyatan {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// lookup of a 32-entry table, j < 32 (no gathers)
template <class D, class VU>
HWY_INLINE hn::VFromD<D> Lookup32(D d, const double *tab, const double *pad, VU j) {
    constexpr size_t T = 32;
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    if (N >= T && N <= 64) {  // table (repeated) in one vector
        return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, hn::BitCast(di, j)));
    }
    if (N >= 2) {  // chunks of 2N entries with TwoTablesLookupLanes, merged by the high index bits
        const hn::RebindToUnsigned<D> du;
        const auto jl = hn::IndicesFromVec(d, hn::BitCast(di, hn::And(j, hn::Set(du, 2 * N - 1))));
        auto r = hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + N), jl);
        for (size_t c = 2 * N; c < T; c += 2 * N) {
            r = hn::IfThenElse(hn::RebindMask(d, hn::Ge(j, hn::Set(du, c))),
                               hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab + c), hn::LoadU(d, tab + c + N), jl), r);
        }
        return r;
    }
    return hn::GatherIndex(d, tab, hn::BitCast(di, j));  // single-lane targets
}

// approximation of 1/y with 16 significant bits: vrcp14pd (native, AVX-512 only) or trunc16(RN(1/y))
template <bool kNative, class D>
HWY_INLINE hn::VFromD<D> RcpApprox(D d, hn::VFromD<D> y) {
#if HWY_TARGET <= HWY_AVX3
    if constexpr (kNative) return hn::ApproximateReciprocal(y);
#endif
    const hn::RebindToUnsigned<D> du;
    const auto r = hn::Div(hn::Set(d, 1.0), y);
    return hn::BitCast(d, hn::And(hn::BitCast(du, r), hn::Set(du, 0xfffffff000000000ULL)));
}

template <bool kNative, class D>
HWY_INLINE hn::VFromD<D> AtanCore(D d, hn::VFromD<D> x) {
    using V = hn::VFromD<D>;
    const hn::RebindToUnsigned<D> du;
    const V one = hn::Set(d, kOne);
    const V ax = hn::Abs(x);
    const auto small = hn::Lt(ax, hn::Set(d, kTabMax));  // a < 7.875

    // table index j and x0 = j/4: a rounded to a multiple of 1/4
    const V S = hn::Set(d, kShifter);
    const V t = hn::Add(ax, S);
    const V x0 = hn::Sub(t, S);
    const auto j = hn::And(hn::BitCast(du, t), hn::Set(du, 31));
    const V diff = hn::IfThenElse(small, hn::Sub(ax, x0), hn::Set(d, kMOne));      // vreducepd, or -1
    const V y = hn::IfThenElse(small, hn::MulAdd(ax, x0, one), hn::Min(hn::Set(d, kYMax), ax));
    const V thi = hn::IfThenElse(small, Lookup32(d, AtanHi(), AtanHiPad(), j), hn::Set(d, kPio2Hi));
    const V tlo = hn::IfThenElse(small, Lookup32(d, AtanLo(), AtanLoPad(), j), hn::Set(d, kPio2Lo));
    const V xx0lo = hn::MulSub(x0, ax, hn::Sub(y, one));  // (x*x0)_low

    // Rcp2 = Rcp (1 + eps)(1 + eps^2), eps = 1 - Y Rcp
    V rcp = RcpApprox<kNative>(d, y);
    const V eps = hn::NegMulAdd(rcp, y, one);
    const V eps2 = hn::Mul(eps, eps);
    rcp = hn::MulAdd(eps, rcp, rcp);
    rcp = hn::MulAdd(eps2, rcp, rcp);

    // R + R_low = Diff/Y
    const V r = hn::Mul(rcp, diff);
    const V eps1 = hn::NegMulAdd(rcp, y, one);
    V rlo = hn::MulAdd(r, eps1, hn::MulSub(rcp, diff, r));
    rlo = hn::IfThenElse(small, hn::NegMulAdd(hn::Mul(xx0lo, rcp), r, rlo), rlo);

    const V r2 = hn::Mul(r, r);
    const V r4 = hn::Mul(r2, r2);
    const V r3 = hn::Mul(r2, r);
    V p45 = hn::MulAdd(hn::Set(d, kA5), r2, hn::Set(d, kA4));
    const V p23 = hn::MulAdd(hn::Set(d, kA3), r2, hn::Set(d, kA2));
    const V p01 = hn::MulAdd(hn::Set(d, kA1), r2, hn::Set(d, kA0));
    p45 = hn::MulAdd(r4, p45, p23);
    V p = hn::MulAdd(r4, p45, p01);
    const V hr = hn::Add(thi, r);
    const V hr_err = hn::Sub(r, hn::Sub(hr, thi));
    p = hn::MulAdd(r3, p, hn::Add(hn::Add(rlo, tlo), hr_err));
    const V res = hn::Xor(hn::Add(p, hr), hn::Xor(x, ax));
    return hn::IfThenElse(hn::IsNaN(x), hn::Add(x, x), res);
}

template <bool kNative>
HWY_INLINE void AtanLoop(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(AtanCore<kNative>(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(AtanCore<kNative>(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
}

void AtanArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AtanLoop<false>(x, y, len); }
void AtanNativeArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AtanLoop<true>(x, y, len); }

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void AtanArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::atan(x[i]);
}
void AtanNativeArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AtanArray(x, y, len); }

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 0: std::atan fallback, 1: SVML algorithm, 2: same, and AtanNativeArray uses vrcp14pd
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return HWY_TARGET <= HWY_AVX3 ? 2 : 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyatan
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyatan {
HWY_EXPORT(AtanArray);
HWY_EXPORT(AtanNativeArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAtan(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AtanArray)(x, y, len); }
static void DispatchNative(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AtanNativeArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyatan

extern "C" {
void hwy_atan_array(const double *x, double *y, size_t len) { hwyatan::DispatchAtan(x, y, len); }
void hwy_atan_native_array(const double *x, double *y, size_t len) { hwyatan::DispatchNative(x, y, len); }
int hwy_atan_variant(void) { return hwyatan::DispatchVariant(); }
size_t hwy_atan_lanes(void) { return hwyatan::DispatchLanes(); }
}
#endif  // HWY_ONCE
