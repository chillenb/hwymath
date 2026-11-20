// Constants for the Highway cosh (hwy_cosh.cc), included once (hwy_cosh.cc itself is re-included per
// target). Everything comes from the SVML data (svml/svml_cosh_d_ha_constants.inc, generated from the
// assembly), with the tables made constexpr so that the interleaved copy below is built at compile time.
#ifndef HWYCOSH_CONSTANTS_H
#define HWYCOSH_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define COSH_TABLE static constexpr
#include "../svml/svml_cosh_d_ha_constants.inc"

namespace hwycosh {

constexpr double kShifter = vShifter, kLog2e = vLog2e, kLn2Hi = vLn2Hi, kLn2Lo = vLn2Lo;
constexpr double kC2 = vC2, kC4 = vC4, kC6 = vC6, kS3 = vS3, kS5 = vS5, kS7 = vS7;
constexpr uint64_t kExpMask = vExpMaskBits, kSignMask = vSignMaskBits;
// callout if the high word of |x| is above vThresholdHi, i.e. bits(|x|) > kThresholdBits (64-bit compare)
constexpr uint64_t kThresholdBits = (uint64_t{vThresholdHi} << 32) | 0xffffffffu;

// the three 16-entry tables (index j = 0..15)
constexpr const double *kTpH = vTpH, *kTpL = vTpL, *kTmH = vTmH;

// The same tables interleaved by index, {TpH, TpL, TmH, 0}[j], for lookups by row (one load per lane and
// a transpose) on targets with 2 or 4 lanes.
struct alignas(64) Rows {
    double v[16 * 4];
    constexpr Rows() : v{} {
        for (size_t j = 0; j < 16; j++) {
            v[4 * j] = vTpH[j];
            v[4 * j + 1] = vTpL[j];
            v[4 * j + 2] = vTmH[j];
            v[4 * j + 3] = 0.0;
        }
    }
};
inline constexpr Rows kRows{};

}  // namespace hwycosh

// the constant macros would clash with identifiers in Highway's headers (e.g. vC2)
#undef COSH_TABLE
#undef vShifter
#undef vC2
#undef vS3
#undef vC4
#undef vS5
#undef vC6
#undef vS7
#undef vLog2e
#undef vLn2Hi
#undef vLn2Lo
#undef vExpMaskBits
#undef vSignMaskBits
#undef vThresholdHi
#undef sTwo
#undef s64OverLn2
#undef sShifter
#undef sP2
#undef sP3
#undef sP4
#undef sP5
#undef sP6
#undef sOverflowX
#undef sLn2o64Hi
#undef sLn2o64Lo
#undef sHuge
#undef sOne
#undef sLargeX

#endif
