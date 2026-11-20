// Constants for the Highway tanh (hwy_tanh.cc), included once (hwy_tanh.cc itself is re-included per
// target). Everything comes from the SVML data (svml/svml_tanh_d_ha_constants.inc, generated from the
// assembly), with the tables made constexpr so that the interleaved copy below is built at compile time.
#ifndef HWYTANH_CONSTANTS_H
#define HWYTANH_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define TANH_TABLE static constexpr
#include "../svml/svml_tanh_d_ha_constants.inc"

namespace hwytanh {

constexpr uint64_t kSignMask = vSignMaskBits, kAbsMask = vAbsMaskBits;
// subinterval index j = clamp(((hi & 0x7ff80000) - 0x3fc00000) >> 19, 0, 15) with hi = bits(x) >> 32
//                     = clamp(u - kIndexBase, 0, 15) with u = (bits(x) >> 51) & 0xfff
constexpr uint32_t kIndexBase = vHiBase >> 19, kIndexMax = (vHiBase + vHiMax) >> 19;
static_assert(vHiMask == 0x7ff80000u && (vHiBase & 0x7ffff) == 0 && vHiMax == 15u << 19, "index layout");
// callout if (hi & 0x7ff80000) > vHugeHi, i.e. bits(|x|) > kThresholdBits (64-bit compare)
constexpr uint64_t kThresholdBits = (uint64_t{vHugeHi | 0x7ffff} << 32) | 0xffffffffu;

constexpr int kTables = 19;  // B, P0 high, P0 low, P1 .. P16
constexpr const double (*kTab)[16] = vTab;

// The same tables interleaved by subinterval, {B, P0h, P0l, P1, ..., P16, 0}[j] (20 doubles per row), for
// lookups by row (one load per lane and group of 4 tables, then a transpose) on targets with 2 or 4 lanes.
constexpr int kRowLen = 20;
struct alignas(64) Rows {
    double v[16 * kRowLen];
    constexpr Rows() : v{} {
        for (int j = 0; j < 16; j++) {
            for (int k = 0; k < kRowLen; k++) v[kRowLen * j + k] = k < kTables ? vTab[k][j] : 0.0;
        }
    }
};
inline constexpr Rows kRows{};

}  // namespace hwytanh

#undef TANH_TABLE
#undef vHiMask
#undef vHiBase
#undef vHiMax
#undef vHugeHi
#undef vSignMaskBits
#undef vAbsMaskBits

#endif
