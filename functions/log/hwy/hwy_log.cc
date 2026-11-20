// Portable (Google Highway) version of SVML's high-accuracy double-precision log (svml/svml_log_d_ha.c, the
// translation of svml_z0_log_d_ha.s), with run-time dispatch to the best target. On every target with true
// FMA (HWY_NATIVE_FMA) it computes the same operations as Intel's AVX-512 code, with the AVX-512-only ones
// emulated exactly (functions/log/hwy/log_hwy-inl.h): vgetexppd/vgetmantpd by integer operations (denormals
// pre-scaled by 2^52), vrcp14pd + vrndscalepd by the exact 16-threshold model of the rounded reciprocal
// (rcp14_thresholds.h), vpermt2pd by TwoTablesLookupLanes (one instruction on AVX-512, two-table lookups
// and selects on narrower vectors). Results are bit-identical to __svml_log8_ha. Lanes with x <= 0,
// +-inf or NaN use the scalar SVML callout (svml_dlog_ha_cout_rare), like the assembly.
// Targets without FMA (SSE2/SSSE3/SSE4): std::log.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dlog_ha_cout_rare(const double *px, double *pr);

#include "functions/log/hwy/hwy_log_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/log/hwy/hwy_log.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
#include "functions/log/hwy/log_hwy-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath_log {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for x <= 0, +-inf, NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void LogFixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(x[k] > 0.0 && x[k] < HUGE_VAL)) svml_dlog_ha_cout_rare(&x[k], &y[k]);
    }
}

// the vector path of __svml_log8_ha (same operations, same order of roundings)
// kDenormals: also correct for denormal x (slower; only used for vectors with special lanes)
template <bool kDenormals, class D, class V = hn::VFromD<D>>
HWY_INLINE V LogCore(D d, V x) {
    using namespace log_e;
    hn::VFromD<hn::RebindToUnsigned<D>> idx;
    V e, m, rcp;
    if (kDenormals) ExpMant(d, x, e, m);
    else ExpMantNormal(d, x, e, m);
    RoundedRcp14(d, m, rcp, idx);
    const auto r = hn::MulSub(rcp, m, hn::Set(d, kOne));  // R = rcp*m - 1
    const auto lt34 = hn::Lt(rcp, hn::Set(d, kThree4));
    const auto p67 = hn::MulAdd(hn::Set(d, kC7), r, hn::Set(d, kC6));
    const auto p89 = hn::MulAdd(hn::Set(d, kC9), r, hn::Set(d, kC8));
    const auto r2 = hn::Mul(r, r);
    const auto p45 = hn::MulAdd(hn::Set(d, kC5), r, hn::Set(d, kC4));
    const auto p69 = hn::MulAdd(r2, p89, p67);
    const auto r4 = hn::Mul(r2, r2);
    const auto p23 = hn::MulAdd(hn::Set(d, kC3), r, hn::Set(d, kC2));
    e = hn::MaskedAddOr(e, lt34, e, hn::Set(d, kOne));  // E + 1 if rcp < 0.75
    const auto th = Lookup16(d, kTh, ThPad(), idx);
    const auto tl = Lookup16(d, kTl, TlPad(), idx);
    const auto p25 = hn::MulAdd(r2, p45, p23);
    const auto high = hn::MulAdd(e, hn::Set(d, kL2H), th);  // E*log(2)_high + Th
    const auto p = hn::MulAdd(r4, p69, p25);
    const auto hr = hn::Add(high, r);
    const auto rh = hn::Sub(hr, high);
    const auto rl = hn::Sub(r, rh);
    const auto q = hn::MulAdd(r2, p, rl);
    const auto lo = hn::MulAdd(hn::Set(d, kL2L), e, tl);  // E*log(2)_low + Tl
    return hn::Add(hn::Add(q, lo), hr);
}

void LogArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto dbl_min = hn::Set(d, 0x1p-1022), inf = hn::Inf(d);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        if (HWY_LIKELY(hn::AllTrue(d, hn::And(hn::Ge(v, dbl_min), hn::Lt(v, inf))))) {
            hn::StoreU(LogCore<false>(d, v), d, y + i);
        } else {  // denormal or special lanes
            hn::StoreU(LogCore<true>(d, v), d, y + i);
            LogFixSpecial(x + i, y + i, N);
        }
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(LogCore<true>(d, v), d, y + i, rem);
        LogFixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm's FMAs would not be exact

void LogArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::log(x[i]);
}

#endif

// Highway's own log (hwy/contrib/math), for comparison
void ContribLogArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(hn::Log(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(hn::Log(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
#else
    for (size_t i = 0; i < len; i++) y[i] = std::log(x[i]);
#endif
}

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm (bit-identical to Intel's AVX-512 code), 0: std::log fallback
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwymath_log
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwymath_log {
HWY_EXPORT(LogArray);
HWY_EXPORT(ContribLogArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchLog(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(LogArray)(x, y, len); }
static void DispatchContrib(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ContribLogArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwymath_log

extern "C" {
void hwy_log_array(const double *x, double *y, size_t len) { hwymath_log::DispatchLog(x, y, len); }
void hwy_contrib_log_array(const double *x, double *y, size_t len) { hwymath_log::DispatchContrib(x, y, len); }
int hwy_log_variant(void) { return hwymath_log::DispatchVariant(); }
size_t hwy_log_lanes(void) { return hwymath_log::DispatchLanes(); }
}
#endif  // HWY_ONCE
