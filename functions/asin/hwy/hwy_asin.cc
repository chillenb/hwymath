// Portable (Google Highway) version of SVML's high-accuracy asin (svml/svml_asin_d_ha.c), with run-time
// dispatch to the best target. Per target:
//   - targets with true FMA (HWY_NATIVE_FMA): the SVML algorithm, operation by operation (all other
//     AVX-512 features have exact equivalents: rounding is round-to-nearest everywhere, vminpd = Min,
//     masked operations = IfThenElse). vrsqrt14pd, whose result bits are implementation-defined, is
//     replaced by trunc16(RN(1/RN(sqrt(Y)))) (16 significant bits, like Intel's, relative error
//     < 2^-15.9 instead of 2^-14), so the results are identical on all these targets and machines
//     (= svml_asin8_ha_portable), and differ from Intel's in the last bit for a small fraction of the
//     inputs with |x| >= 0.5 (see README.md);
//     AsinNativeArray uses vrsqrt14pd itself on AVX-512 targets (bit-identical to Intel's code on the
//     same CPU) and is the same as AsinArray elsewhere;
//   - targets without FMA (SSE2/SSSE3/SSE4): std::asin.
// Lanes with |x| > 1 (and +-inf) use the scalar SVML callout (svml_dasin_ha_cout_rare); NaN lanes return
// the quieted input, like the SVML vector path.
// Based on Intel SVML (BSD-3-Clause, Copyright (C) 2023 Intel Corporation); same license.

#include <cmath>
#include <cstddef>
#include <cstdint>

extern "C" int svml_dasin_ha_cout_rare(const double *px, double *pr);

#include "functions/asin/hwy/hwy_asin_constants.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/asin/hwy/hwy_asin.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace hwyasin {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64

// scalar fix-up for lanes with |x| > 1 (rare: kept out of line)
static HWY_NOINLINE __attribute__((cold)) void FixSpecial(const double *x, double *y, size_t n) {
    for (size_t k = 0; k < n; k++) {
        if (std::fabs(x[k]) > kOne) svml_dasin_ha_cout_rare(&x[k], &y[k]);
    }
}

// approximation of 1/sqrt(y) with 16 significant bits: vrsqrt14pd (native, AVX-512 only) or the
// deterministic trunc16(RN(1/RN(sqrt(y))))
template <bool kNative, class D>
HWY_INLINE hn::VFromD<D> RsqrtApprox(D d, hn::VFromD<D> y) {
#if HWY_TARGET <= HWY_AVX3
    if constexpr (kNative) return hn::ApproximateReciprocalSqrt(y);
#endif
    const hn::RebindToUnsigned<D> du;
    const auto r = hn::Div(hn::Set(d, 1.0), hn::Sqrt(y));
    return hn::BitCast(d, hn::And(hn::BitCast(du, r), hn::Set(du, 0xfffffff000000000ULL)));
}

template <bool kNative, class D>
HWY_INLINE hn::VFromD<D> AsinCore(D d, hn::VFromD<D> x) {
    using V = hn::VFromD<D>;
    const V half = hn::Set(d, kHalf);
    const V ax = hn::Abs(x);
    const V x2 = hn::Mul(ax, ax);
    const auto big = hn::Not(hn::Lt(ax, half));  // a >= 0.5

    const V y = hn::NegMulAdd(half, ax, half);  // Y = 0.5 - 0.5 a
    V r = RsqrtApprox<kNative>(d, y);
    r = hn::IfThenZeroElse(hn::Lt(y, hn::Set(d, kSmallY)), r);
    const V rr = hn::Min(x2, y);  // R

    // 2 sqrt(Y) = Sh - Sl
    const V y2 = hn::Add(y, y);
    const V r2 = hn::Mul(r, r);
    const V sh = hn::Mul(y2, r);
    const V e = hn::MulSub(r2, y2, hn::Set(d, kTwo));
    const V sherr = hn::MulSub(r, y2, sh);
    V sp = hn::MulAdd(hn::Set(d, kSqrtC3), e, hn::Set(d, kSqrtC2));
    sp = hn::MulAdd(e, sp, hn::Set(d, kSqrtC1));
    sp = hn::MulAdd(e, sp, hn::Set(d, kSqrtC0));
    const V sl = hn::MulSub(hn::Mul(sh, e), sp, sherr);
    const V s = hn::IfThenElse(big, hn::Sub(sl, sh), ax);  // -2 sqrt(Y), or a

    // polynomial
    const V rr2 = hn::Mul(rr, rr);
    const V rr4 = hn::Mul(rr2, rr2);
    const V q910 = hn::MulAdd(hn::Set(d, kC10), rr, hn::Set(d, kC9));
    const V q34 = hn::MulAdd(hn::Set(d, kC4), rr, hn::Set(d, kC3));
    const V q56 = hn::MulAdd(hn::Set(d, kC6), rr, hn::Set(d, kC5));
    const V q1112 = hn::MulAdd(rr2, hn::MulAdd(hn::Set(d, kC12), rr, hn::Set(d, kC11)), q910);
    const V q78 = hn::MulAdd(rr2, hn::MulAdd(hn::Set(d, kC8), rr, hn::Set(d, kC7)), q56);
    V q = hn::MulAdd(rr4, q1112, q78);
    q = hn::MulAdd(rr2, q, q34);
    q = hn::MulAdd(rr, q, hn::Set(d, kC2));
    q = hn::MulAdd(rr, q, hn::Set(d, kC1));
    V pr = hn::Mul(q, rr);

    // pi/2 - Sh with its rounding error
    const V pio2 = hn::Set(d, kPio2Hi);
    const V hi = hn::Sub(pio2, sh);
    const V lo = hn::Add(hn::Set(d, kPio2Lo), sl);
    const V err = hn::Sub(sh, hn::Sub(pio2, hi));
    const V low = hn::IfThenElseZero(big, hn::Sub(lo, err));
    pr = hn::MulAdd(s, pr, low);
    const V res = hn::Add(hn::IfThenElse(big, hi, s), pr);
    const V signed_res = hn::Xor(res, hn::Xor(x, ax));
    return hn::IfThenElse(hn::IsNaN(x), hn::Add(x, x), signed_res);
}

template <bool kNative>
HWY_INLINE void AsinLoop(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    const auto one = hn::Set(d, kOne);
    size_t i = 0;
    for (; i + N <= len; i += N) {
        const auto v = hn::LoadU(d, x + i);
        hn::StoreU(AsinCore<kNative>(d, v), d, y + i);
        if (HWY_UNLIKELY(!hn::AllFalse(d, hn::Lt(one, hn::Abs(v))))) FixSpecial(x + i, y + i, N);
    }
    if (i < len) {
        const size_t rem = len - i;
        const auto v = hn::LoadN(d, x + i, rem);
        hn::StoreN(AsinCore<kNative>(d, v), d, y + i, rem);
        FixSpecial(x + i, y + i, rem);
    }
}

void AsinArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AsinLoop<false>(x, y, len); }
void AsinNativeArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AsinLoop<true>(x, y, len); }

#else  // no fused multiply-add: the algorithm's error-free steps would not be exact

void AsinArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
    for (size_t i = 0; i < len; i++) y[i] = std::asin(x[i]);
}
void AsinNativeArray(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) { AsinArray(x, y, len); }

#endif

size_t LanesF64() { return hn::Lanes(hn::ScalableTag<double>()); }

int Variant() {  // 0: std::asin fallback, 1: SVML algorithm, 2: same, and AsinNativeArray uses vrsqrt14pd
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    return HWY_TARGET <= HWY_AVX3 ? 2 : 1;
#else
    return 0;
#endif
}

}  // namespace HWY_NAMESPACE
}  // namespace hwyasin
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyasin {
HWY_EXPORT(AsinArray);
HWY_EXPORT(AsinNativeArray);
HWY_EXPORT(Variant);
HWY_EXPORT(LanesF64);
// dispatch inside the namespace: with a single (static) target, HWY_DYNAMIC_DISPATCH(F) becomes N_<TARGET>::F
static void DispatchAsin(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AsinArray)(x, y, len); }
static void DispatchNative(const double *x, double *y, size_t len) { HWY_DYNAMIC_DISPATCH(AsinNativeArray)(x, y, len); }
static int DispatchVariant() { return HWY_DYNAMIC_DISPATCH(Variant)(); }
static size_t DispatchLanes() { return HWY_DYNAMIC_DISPATCH(LanesF64)(); }
}  // namespace hwyasin

extern "C" {
void hwy_asin_array(const double *x, double *y, size_t len) { hwyasin::DispatchAsin(x, y, len); }
void hwy_asin_native_array(const double *x, double *y, size_t len) { hwyasin::DispatchNative(x, y, len); }
int hwy_asin_variant(void) { return hwyasin::DispatchVariant(); }
size_t hwy_asin_lanes(void) { return hwyasin::DispatchLanes(); }
}
#endif  // HWY_ONCE
