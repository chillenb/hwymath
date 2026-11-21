// Portable (Google Highway) version of SVML's high-accuracy double-precision log1p (svml/svml_log1p_d_ha.c, the
// translation of svml_z0_log1p_d_ha.s), with run-time dispatch to the best target. On every target with true
// FMA (HWY_NATIVE_FMA) it computes the same operations as Intel's AVX-512 code, with the AVX-512-only ones
// emulated exactly (functions/log/hwy/log_hwy-inl.h): vgetexppd/vgetmantpd by integer operations (1 + x is normal),
// vscalefpd by two exact power-of-2 multiplications, vrcp14pd + vrndscalepd by the exact 16-threshold model of the rounded reciprocal
// (rcp14_thresholds.h), vpermt2pd by TwoTablesLookupLanes (one instruction on AVX-512, two-table lookups
// and selects on narrower vectors). Results are bit-identical to __svml_log1p8_ha. Lanes with 1 + x <= 0,
// +-inf or NaN use the scalar SVML callout (svml_dlog1p_ha_cout_rare), like the assembly.
// On AVX-512 targets (HWYMATH_NATIVE_AVX512) the instructions themselves replace these emulations.
// Targets without FMA (SSE2/SSSE3/SSE4): std::log.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dlog1p_ha_cout_rare(const double *px, double *pr);

#include "functions/log1p/hwy/hwy_log1p_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/log1p/hwy/hwy_log1p.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwy/contrib/math/math-inl.h"
#include "functions/log/hwy/log_hwy-inl.h"
#include "hwymath/hwy/avx512_native-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwymath_log {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for 1 + x <= 0, +-inf, NaN (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void Log1pFixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (!(x[k] + 1.0 > 0.0 && x[k] + 1.0 < HUGE_VAL)) svml_dlog1p_ha_cout_rare(&x[k], &y[k]);
    }
}

// the vector path of __svml_log1p8_ha (same operations, same order of roundings)
template <class D, class V = hn::VFromD<D>>
HWY_INLINE V Log1pCore(D d, V x) {
    using namespace log1p_e;
    const hn::RebindToUnsigned<D> du;
    hn::VFromD<hn::RebindToUnsigned<D>> idx;
    V e, m, rcp;
    // 1 + x = Xh + Xl
    const auto one = hn::Set(d, kOne);
    const auto mx = hn::Max(one, x);
    const auto mn = hn::Min(one, x);
    const auto xh = hn::Add(one, x);
    const auto sgn = hn::And(hn::BitCast(du, x), hn::Set(du, kSignMask));
    ExpMantNormal(d, xh, e, m);
    const auto xl = hn::Sub(mn, hn::Sub(xh, mx));
    RoundedRcp14(d, m, rcp, idx);
#if HWYMATH_NATIVE_AVX512
    const auto xls = hwymath::HWY_NAMESPACE::Scalef(xl, hn::Neg(e));  // Xl*2^-E: vxorpd, vscalefpd
#else
    // Xl*2^-E (vscalefpd): exact multiplication by 2^(60-E), then one rounding in the multiplication by 2^-60
    const auto p60 = hn::ShiftLeft<52>(hn::Sub(hn::Set(du, 2106), hn::ShiftRight<52>(hn::BitCast(du, xh))));
    const auto xls = hn::Mul(hn::Mul(xl, hn::BitCast(d, p60)), hn::Set(d, 0x1p-60));
#endif
    const auto lt34 = hn::Lt(rcp, hn::Set(d, kThree4));
    const auto r = hn::MulSub(rcp, m, one);  // R = rcp*m - 1
    const auto rxl = hn::Mul(rcp, xls);
    e = hn::MaskedAddOr(e, lt34, e, one);  // E + 1 if rcp < 0.75
    const auto rr = hn::Add(r, rxl);       // R'
    const auto r2 = hn::Mul(rr, rr);
    const auto p89 = hn::MulAdd(hn::Set(d, kC9), rr, hn::Set(d, kC8));
    const auto rl = hn::Sub(rxl, hn::Sub(rr, r));
    const auto tl = Lookup16(d, kTl, TlPad(), idx);
    const auto th = Lookup16(d, kTh, ThPad(), idx);
    const auto high = hn::MulAdd(hn::Set(d, kL2H), e, th);  // E*log(2)_high + Th
    const auto tlrl = hn::Add(tl, rl);
    const auto p67 = hn::MulAdd(hn::Set(d, kC7), rr, hn::Set(d, kC6));
    const auto p23 = hn::MulAdd(hn::Set(d, kC3), rr, hn::Set(d, kC2));
    const auto hr = hn::Add(high, rr);
    const auto p69 = hn::MulAdd(r2, p89, p67);
    const auto p45 = hn::MulAdd(hn::Set(d, kC5), rr, hn::Set(d, kC4));
    const auto r4 = hn::Mul(r2, r2);
    const auto p25 = hn::MulAdd(r2, p45, p23);
    const auto rh = hn::Sub(hr, high);
    const auto p = hn::MulAdd(r4, p69, p25);
    const auto rl2 = hn::Sub(rr, rh);
    const auto q = hn::MulAdd(r2, p, rl2);
    const auto lo = hn::MulAdd(hn::Set(d, kL2L), e, tlrl);  // E*log(2)_low + Tl + Rl
    const auto res = hn::Add(hr, hn::Add(q, lo));
    return hn::BitCast(d, hn::Or(hn::BitCast(du, res), sgn));  // log1p(-0) = -0
}

template <class D, class V = hn::VFromD<D>>
HWY_INLINE auto Log1pSpecial(D d, V x) {  // 1 + x <= 0, +-inf, NaN
#if HWYMATH_NATIVE_AVX512
    // vfpclasspd $0xdf on 1 + x: every class but positive normal and denormal (1 + x is never denormal)
    return hwymath::HWY_NAMESPACE::FpClass<0xdf>(hn::Add(x, hn::Set(d, 1.0)));
#else
    return hn::Not(hn::And(hn::Gt(hn::Add(x, hn::Set(d, 1.0)), hn::Zero(d)), hn::Lt(hn::Add(x, hn::Set(d, 1.0)), hn::Inf(d))));
#endif
}

void Log1pArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(Log1pCore(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllFalse(d, Log1pSpecial(d, v)))) Log1pFixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(Log1pCore(d, v), d, y + i, rem);
        Log1pFixSpecial(x + i, y + i, rem);
    }
}

#else  // no fused multiply-add: the algorithm's FMAs would not be exact

void Log1pArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::log1p(x[i]);
}

#endif

// Highway's own log1p (hwy/contrib/math), for comparison
void ContribLog1pArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) hn::StoreU(hn::Log1p(d, hn::LoadU(d, x + i)), d, y + i);
    if (i < len) hn::StoreN(hn::Log1p(d, hn::LoadN(d, x + i, len - i)), d, y + i, len - i);
#else
    for (size_t i = 0; i < len; i++) y[i] = std::log1p(x[i]);
#endif
}

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 1: SVML algorithm (bit-identical to Intel's AVX-512 code), 0: std::log1p fallback
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
HWY_EXPORT(Log1pArray);
HWY_EXPORT(ContribLog1pArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchLog1p(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(Log1pArray)(x, y, len); }
static void DispatchContrib(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(ContribLog1pArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwymath_log

extern "C" {
void hwy_log1p_array(const double *x, double *y, size_t len) { hwymath_log::DispatchLog1p(x, y, len); }
void hwy_contrib_log1p_array(const double *x, double *y, size_t len) { hwymath_log::DispatchContrib(x, y, len); }
int hwy_log1p_variant(void) { return hwymath_log::DispatchVariant(); }
size_t hwy_log1p_lanes(void) { return hwymath_log::DispatchLanes(); }
}
#endif  // HWY_ONCE
