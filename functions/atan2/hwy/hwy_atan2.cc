// Portable (Google Highway) version of SVML's high-accuracy atan2 (svml/svml_atan2_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation by operation (all other
//     AVX-512 features have exact equivalents: rounding is round-to-nearest everywhere, masked
//     operations = IfThenElse, the 32-bit high-word range check = two compares of |v|);
//     vrcp14pd is emulated bit-exactly (common/include/hwymath/hwy/avx512_emul-inl.h, after Intel's reference RECIP14.c) on
//     targets without AVX-512 and is the instruction itself on AVX-512 targets, so the results are
//     bit-identical to Intel's AVX-512 code (svml_atan28_ha) on all these targets;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::atan2.
// Lanes where |x| or |y| is outside [2^-1020, 2^993) (zeros, denormals, huge values, inf, NaN) are done in
// scalar code: zero arguments (no NaN) like the vector code's zero branch (svml_atan2_ha_zero_args), the
// others by the scalar SVML callout (svml_datan2_ha_cout_rare).
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_datan2_ha_cout_rare(const double *py, const double *px, double *pr);
extern "C" double svml_atan2_ha_zero_args(double y, double x);

#include "functions/atan2/hwy/hwy_atan2_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/atan2/hwy/hwy_atan2.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "hwymath/hwy/avx512_emul-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyatan2 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

static inline bool InRange(double v) {  // false for NaN
    const double a = std::fabs(v);
    return a >= kMinArg && a < kMaxArg;
}

// scalar fix-up for lanes outside the vector path's range (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *y, const double *x, double *r, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (InRange(x[k]) && InRange(y[k])) continue;
        if ((x[k] == 0.0 || y[k] == 0.0) && !std::isnan(x[k]) && !std::isnan(y[k])) {
            r[k] = svml_atan2_ha_zero_args(y[k], x[k]);
        } else {
            svml_datan2_ha_cout_rare(&y[k], &x[k], &r[k]);
        }
    }
}

// vrcp14pd: the instruction on AVX-512 targets, elsewhere its bit-exact emulation (common/include/hwymath/hwy/avx512_emul-inl.h; its
// domain covers every lane whose result is used)
template <class D>
HWY_INLINE hn::VFromD<D> Rcp14(D d, hn::VFromD<D> b) {
    return hwymath::HWY_NAMESPACE::Rcp14(d, b);
}

template <class D>
HWY_INLINE hn::VFromD<D> Atan2Core(D d, hn::VFromD<D> y, hn::VFromD<D> x) {
    using V = hn::VFromD<D>;
    const V one = hn::Set(d, kOne);
    const V ax = hn::Abs(x), ay = hn::Abs(y);
    const auto xneg = hn::Lt(x, hn::Zero(d));
    const V sgnx = hn::Xor(x, ax), sgny = hn::Xor(y, ay);

    // interval of |y/x|
    const auto k5 = hn::Lt(hn::Mul(ax, hn::Set(d, kT7)), ay);
    const auto k1 = hn::Lt(hn::Mul(ax, hn::Set(d, kT11)), ay);
    const auto k2 = hn::Lt(hn::Mul(ax, hn::Set(d, kT19)), ay);
    const auto k3 = hn::Lt(hn::Mul(ax, hn::Set(d, kT39)), ay);
    const V dd = hn::IfThenElse(k2, hn::IfThenElse(k3, one, hn::Set(d, kOneHalf)),
                                hn::IfThenElse(k1, one, hn::Set(d, kHalf)));
    const V hi = hn::IfThenElse(k2, hn::IfThenElse(k3, hn::Set(d, kPio2Hi), hn::Set(d, kAtan15Hi)),
                                hn::IfThenElse(k1, hn::Set(d, kPio4Hi), hn::Set(d, kAtan05Hi)));
    const V lo = hn::IfThenElse(k2, hn::IfThenElse(k3, hn::Set(d, kPio2Lo), hn::Set(d, kAtan15Lo)),
                                hn::IfThenElse(k1, hn::Set(d, kPio4Lo), hn::Set(d, kAtan05Lo)));
    V b = hn::IfThenZeroElse(k3, ax);
    b = hn::IfThenElse(k5, hn::MulAdd(dd, ay, b), b);
    V a = hn::IfThenZeroElse(k3, ay);
    a = hn::IfThenElse(k5, hn::NegMulAdd(dd, ax, a), a);

    // 1/b ~ r2, by two refinements of r0; R + Rlo = a/b
    V r = Rcp14(d, b);
    r = hn::MulAdd(hn::NegMulAdd(b, r, one), r, r);
    r = hn::MulAdd(hn::NegMulAdd(b, r, one), r, r);
    const V rr = hn::Mul(r, a);
    const V r2 = hn::Mul(rr, rr);
    const V r4 = hn::Mul(r2, r2);
    V rlo = hn::Mul(hn::NegMulAdd(rr, b, a), r);
    rlo = hn::IfThenElse(k5, hn::Add(rlo, lo), rlo);

    // polynomial
    V pa = hn::MulAdd(hn::Set(d, kC12), r4, hn::Set(d, kC10));
    V pb = hn::MulAdd(hn::Set(d, kC11), r4, hn::Set(d, kC9));
    pa = hn::MulAdd(r4, pa, hn::Set(d, kC8));
    pb = hn::MulAdd(r4, pb, hn::Set(d, kC7));
    pa = hn::MulAdd(r4, pa, hn::Set(d, kC6));
    pb = hn::MulAdd(r4, pb, hn::Set(d, kC5));
    pa = hn::MulAdd(r4, pa, hn::Set(d, kC4));
    pb = hn::MulAdd(r4, pb, hn::Set(d, kC3));
    pa = hn::MulAdd(r4, pa, hn::Set(d, kC2));
    pb = hn::MulAdd(r4, pb, hn::Set(d, kC1));
    const V p = hn::Mul(hn::MulAdd(r2, pa, pb), r2);

    // low parts, then R, then atan(c)_hi; sign of x, + pi for x < 0, sign of y
    rlo = hn::IfThenElse(xneg, hn::Add(rlo, hn::Xor(sgnx, hn::Set(d, kPiLo))), rlo);
    V res = hn::Add(rr, hn::MulAdd(rr, p, rlo));
    res = hn::IfThenElse(k5, hn::Add(res, hi), res);
    res = hn::Xor(res, sgnx);
    res = hn::IfThenElse(xneg, hn::Add(res, hn::Set(d, kPiHi)), res);
    return hn::Or(res, sgny);
}

template <class D>
HWY_INLINE bool AnySpecial(D d, hn::VFromD<D> y, hn::VFromD<D> x) {
    const auto lo = hn::Set(d, kMinArg), hi = hn::Set(d, kMaxArg);
    const auto ax = hn::Abs(x), ay = hn::Abs(y);
    const auto ok = hn::And(hn::And(hn::Ge(ax, lo), hn::Lt(ax, hi)), hn::And(hn::Ge(ay, lo), hn::Lt(ay, hi)));
    return !hn::AllTrue(d, ok);
}

HWY_INLINE void Atan2Loop(const double *HWY_RESTRICT y, const double *HWY_RESTRICT x, double *HWY_RESTRICT r,
                          size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto vy = hn::LoadU(d, y + i), vx = hn::LoadU(d, x + i);
        hn::StoreU(Atan2Core(d, vy, vx), d, r + i);
        if (HWY_UNLIKELY(AnySpecial(d, vy, vx))) FixSpecial(y + i, x + i, r + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto vy = hn::LoadN(d, y + i, rem), vx = hn::LoadN(d, x + i, rem);
        hn::StoreN(Atan2Core(d, vy, vx), d, r + i, rem);
        FixSpecial(y + i, x + i, r + i, rem);
    }
}

void Atan2Array(const double *HWY_RESTRICT y, const double *HWY_RESTRICT x, double *HWY_RESTRICT r, size_t len) {
    Atan2Loop(y, x, r, len);
}

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void Atan2Array(const double *HWY_RESTRICT y, const double *HWY_RESTRICT x, double *HWY_RESTRICT r, size_t len) {
    for (size_t i = 0; i < len; i++) r[i] = std::atan2(y[i], x[i]);
}

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 0: std fallback, 1: SVML algorithm (bit-identical to Intel's)
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyatan2
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyatan2 {
HWY_EXPORT(Atan2Array);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAtan2(const double *y, const double *x, double *r, size_t len) {
    HWY_DYNAMIC_DISPATCH(Atan2Array)(y, x, r, len);
}
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyatan2

extern "C" {
void hwy_atan2_array(const double *y, const double *x, double *r, size_t len) { hwyatan2::DispatchAtan2(y, x, r, len); }
int hwy_atan2_variant(void) { return hwyatan2::DispatchVariant(); }
size_t hwy_atan2_lanes(void) { return hwyatan2::DispatchLanes(); }
}
#endif  // HWY_ONCE
