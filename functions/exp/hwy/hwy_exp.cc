// Portable (Google Highway) version of the extra-precision exp (avx2/exp_d_avx2_hp.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the hp algorithm; table size 2^K with K = 3 on AVX-512
//     (8 entries = one vector) and K = 2 elsewhere (4 entries: one vector on AVX2, two on 128-bit
//     targets such as NEON/SSE-class; the lookup adapts to the vector length at run time, also for
//     SVE/RVV);
//   - targets without FMA (SSE2/SSSE3/SSE4): the error-free steps would not be exact, so std::exp.
// Lanes with |x| >= 707.703... or NaN use the scalar SVML callout (svml_dexp_ha_cout_rare).
// Also exports Highway's own hwy/contrib/math Exp for comparison.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dexp_ha_cout_rare(const double *px, double *pr);

#include "functions/exp/hwy/hwy_exp_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/exp/hwy/hwy_exp.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyexp {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for lanes the vector path does not handle: |x| >= threshold or NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(std::fabs(x[k]) < kThreshold)) svml_dexp_ha_cout_rare(&x[k], &y[k]);
    }
}

#if HWY_TARGET <= HWY_AVX3
constexpr int kK = 3;
#else
constexpr int kK = 2;
#endif
using C = ExpConsts<kK>;

// lookup of a 2^kK-entry table; j < 2^kK
template <class D, class VU>
HWY_INLINE hn::VFromD<D> Lookup(D d, const double *tab, const double *pad, VU j) {
    constexpr size_t T = size_t{1} << kK;
    const hn::RebindToSigned<D> di;
    const size_t N = hn::Lanes(d);
    if (N >= T && N <= 64) {  // table (repeated) in one vector
        return hn::TableLookupLanes(hn::LoadU(d, pad), hn::IndicesFromVec(d, hn::BitCast(di, j)));
    }
    if (2 * N == T) {  // two vectors
        return hn::TwoTablesLookupLanes(d, hn::LoadU(d, tab), hn::LoadU(d, tab + N),
                                        hn::IndicesFromVec(d, hn::BitCast(di, j)));
    }
    return hn::GatherIndex(d, tab, hn::BitCast(di, j));  // e.g. single-lane targets
}

template <class D>
HWY_INLINE hn::VFromD<D> ExpCore(D d, hn::VFromD<D> x) {
    const hn::RebindToUnsigned<D> du;
    const auto S = hn::Set(d, C::kShifter);
    const auto t = hn::MulAdd(x, hn::Set(d, kLog2e), S);  // RN(x*log2(e)) to a multiple of 2^-K
    const auto tb = hn::BitCast(du, t);
    const auto n = hn::Sub(t, S);
    const auto j = hn::And(tb, hn::Set(du, (uint64_t{1} << kK) - 1));
    const auto th = Lookup(d, C::Th(), C::ThPad(), j);
    const auto tlr = Lookup(d, C::Tlr(), C::TlrPad(), j);

    const auto rhi = hn::NegMulAdd(n, hn::Set(d, C::kLn2Hi), x);  // exact
    const auto rlo = hn::Mul(n, hn::Set(d, -C::kLn2Lo));
    const auto r = hn::Add(rhi, rlo);
    const auto rc = hn::Add(hn::Sub(rhi, r), rlo);

    const auto r2 = hn::Mul(r, r);
    const auto q23 = hn::MulAdd(hn::Set(d, C::C3), r, hn::Set(d, C::C2));
    const auto q45 = hn::MulAdd(hn::Set(d, C::C5), r, hn::Set(d, C::C4));
    const auto q67 = hn::MulAdd(hn::Set(d, C::C7), r, hn::Set(d, C::C6));
    hn::VFromD<D> q6;
    if constexpr (C::kDeg == 9) {
        q6 = hn::MulAdd(hn::MulAdd(hn::Set(d, C::C9), r, hn::Set(d, C::C8)), r2, q67);
    } else {
        q6 = hn::MulAdd(hn::Set(d, C::C8), r2, q67);
    }
    const auto q = hn::MulAdd(hn::MulAdd(q6, r2, q45), r2, q23);

    const auto small = hn::MulAdd(r2, q, hn::Add(rc, hn::MulAdd(tlr, r, tlr)));
    const auto s = hn::MulAdd(th, r, th);
    const auto e = hn::MulSub(th, r, hn::Sub(s, th));
    const auto y = hn::Add(s, hn::MulAdd(th, small, e));

    const auto m = hn::ShiftRight<kK>(hn::Sub(tb, hn::Set(du, C::kShifterBase)));
    return hn::Mul(y, hn::BitCast(d, hn::ShiftLeft<52>(m)));
}

void ExpArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto thr = hn::Set(d, kThreshold);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(ExpCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllFalse(d, hn::Not(hn::Lt(hn::Abs(v), thr))))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(ExpCore(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the error-free transformations would not be exact

void ExpArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::exp(x[i]);
}

#endif

// Highway's own exp (hwy/contrib/math), for comparison
void ContribExpArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(hn::Exp(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(hn::Exp(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
#else
    for (size_t i = 0; i < len; i++) y[i] = std::exp(x[i]);
#endif
}

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 3/2: hp algorithm with 2^K table, 0: std::exp fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return kK;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyexp
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyexp {
HWY_EXPORT(ExpArray);
HWY_EXPORT(ContribExpArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchExp(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ExpArray)(x, y, len); }
static void DispatchContrib(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ContribExpArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyexp

extern "C" {
void hwy_exp_array(const double *x, double *y, size_t len) { hwyexp::DispatchExp(x, y, len); }
void hwy_contrib_exp_array(const double *x, double *y, size_t len) { hwyexp::DispatchContrib(x, y, len); }
int hwy_exp_variant(void) { return hwyexp::DispatchVariant(); }
size_t hwy_exp_lanes(void) { return hwyexp::DispatchLanes(); }
}
#endif  // HWY_ONCE
