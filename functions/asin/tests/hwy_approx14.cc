// Array versions of the Highway emulations in hwy/approx14-inl.h, for test_approx14 (after
// functions/asin/tests/hwy_approx14.cc, branch port/invhyperbolic).
#include <cstddef>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "functions/asin/tests/hwy_approx14.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"
#include "functions/asin/hwy/approx14-inl.h"

HWY_BEFORE_NAMESPACE();
namespace hwyinvtrig_approx14 {
namespace HWY_NAMESPACE {
namespace hn = hwy::HWY_NAMESPACE;

template <int kOp>
void Apply(const double *HWY_RESTRICT x, double *HWY_RESTRICT y, size_t len) {
#if HWY_NATIVE_FMA && HWY_HAVE_FLOAT64
    const hn::ScalableTag<double> d;
    const size_t N = hn::Lanes(d);
    for (size_t i = 0; i < len; i += N) {
        const size_t n = len - i < N ? len - i : N;
        const auto v = hn::LoadN(d, x + i, n);
        hn::VFromD<decltype(d)> r;
        if constexpr (kOp == 0) r = hwyinvtrig::HWY_NAMESPACE::Rcp14(d, v);
        else r = hwyinvtrig::HWY_NAMESPACE::Rsqrt14(d, v);
        hn::StoreN(r, d, y + i, n);
    }
#else
    (void)x, (void)y, (void)len;
#endif
}
void Rcp14Array(const double *x, double *y, size_t n) { Apply<0>(x, y, n); }
void Rsqrt14Array(const double *x, double *y, size_t n) { Apply<1>(x, y, n); }
int Emulated() { return HWY_NATIVE_FMA && HWY_HAVE_FLOAT64; }

}  // namespace HWY_NAMESPACE
}  // namespace hwyinvtrig_approx14
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace hwyinvtrig_approx14 {
HWY_EXPORT(Rcp14Array);
HWY_EXPORT(Rsqrt14Array);
HWY_EXPORT(Emulated);
static void D0(const double *x, double *y, size_t n) { HWY_DYNAMIC_DISPATCH(Rcp14Array)(x, y, n); }
static void D1(const double *x, double *y, size_t n) { HWY_DYNAMIC_DISPATCH(Rsqrt14Array)(x, y, n); }
static int D3() { return HWY_DYNAMIC_DISPATCH(Emulated)(); }
}  // namespace hwyinvtrig_approx14
extern "C" {
void hwyinvtrig_rcp14_array(const double *x, double *y, size_t n) { hwyinvtrig_approx14::D0(x, y, n); }
void hwyinvtrig_rsqrt14_array(const double *x, double *y, size_t n) { hwyinvtrig_approx14::D1(x, y, n); }
int hwyinvtrig_approx14_emulated(void) { return hwyinvtrig_approx14::D3(); }
}
#endif
