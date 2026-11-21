// The rounded reciprocal of SVML's log family: rcp = roundscale(vrcp14pd(m), 2^-5, RN) = k/32 for a mantissa m
// in [1, 2). vrcp14pd's low bits are implementation-defined, but after the rounding, k depends on m only
// through these 16 thresholds (all multiples of 2^-16):
//     k(m) = j      for kRcpThreshold[j - 16] <= m < kRcpThreshold[j - 15]   (j = 16..31)
//     k(m) = 32     for m < kRcpThreshold[15]
// i.e. k = 32 - #{thresholds <= m}. kRcpThreshold[j - 16] is the smallest m with k(m) <= j. The exact
// mathematical boundaries would be 64/(2j+1); vrcp14pd's error (< 2^-14) moves them by up to ~2 units of 2^-16.
// Found by bisection on Intel SDE 10.13 (-spr); functions/log/tools/verify_rcp14.c proves them from Intel's
// reference algorithm (common/include/hwymath/rcp14_ref.h, all 2^16 fraction prefixes) and checks them against vrcp14pd.
#ifndef HWYMATH_LOG_RCP14_THRESHOLDS_H
#define HWYMATH_LOG_RCP14_THRESHOLDS_H

#include <cstdint>

namespace hwymath_log {

constexpr double kRcpThreshold[16] = {
    0x1.f08p+0,    // k = 16 above, 17 below   (64/33 = 0x1.f07c1f...p+0)
    0x1.d41cp+0,   // k = 17 / 18              (64/35)
    0x1.badp+0,    // k = 18 / 19              (64/37)
    0x1.a41ap+0,   // k = 19 / 20              (64/39)
    0x1.8f9cp+0,   // k = 20 / 21              (64/41)
    0x1.7d07p+0,   // k = 21 / 22              (64/43)
    0x1.6c15p+0,   // k = 22 / 23              (64/45)
    0x1.5c99p+0,   // k = 23 / 24              (64/47)
    0x1.4e5fp+0,   // k = 24 / 25              (64/49)
    0x1.4143p+0,   // k = 25 / 26              (64/51)
    0x1.3523p+0,   // k = 26 / 27              (64/53)
    0x1.29e6p+0,   // k = 27 / 28              (64/55)
    0x1.1f7p+0,    // k = 28 / 29              (64/57)
    0x1.15b4p+0,   // k = 29 / 30              (64/59)
    0x1.0c97p+0,   // k = 30 / 31              (64/61)
    0x1.040fp+0,   // k = 31 above, 32 below   (64/63)
};

// Emulation (log_hwy-inl.h: RoundedRcp14), with integer operations on the high 32 bits h of m:
// every bin of m with the same top 5 mantissa bits b (width 1/32) contains at most one threshold (they are
// at least 0.0333 apart), so with a 32-entry table kRcpBins[b] = hi32(threshold in the bin) | (k above it - 16)
// (no threshold: hi32(1.0) | (k of the bin - 16)),
//     k - 16 = (kRcpBins[b] & 15) + ((h | 15) < kRcpBins[b])
// (thresholds are multiples of 2^-16, so their high words have 4 zero low bits). This depends on h only:
// verify_rcp14 checks it for all 2^20 values of h.
constexpr uint32_t Hi32(double t) { return 0x3ff00000u | (uint32_t)((t - 1.0) * 65536.0) << 4; }
constexpr int StepK(double m) {  // k(m) of the threshold model
    int k = 32;
    for (int j = 15; j >= 0; j--) if (m >= kRcpThreshold[j]) k = 16 + j;
    return k;
}
struct RcpBins {
    uint32_t v[32];
    constexpr RcpBins() : v{} {
        for (int b = 0; b < 32; b++) {
            const double lo = 1.0 + b / 32.0, hi = 1.0 + (b + 1) / 32.0;
            v[b] = 0x3ff00000u | (uint32_t)(StepK(lo) - 16);
            for (int j = 0; j < 16; j++) {
                if (kRcpThreshold[j] > lo && kRcpThreshold[j] < hi) v[b] = Hi32(kRcpThreshold[j]) | (uint32_t)j;
            }
        }
    }
};
constexpr RcpBins kRcpBins{};

}  // namespace hwymath_log

#endif
