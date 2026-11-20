// Constants for the Highway pow (hwy_pow.cc), included once (hwy_pow.cc itself is re-included per target):
// the SVML tables and constants as constexpr (from ../svml/svml_pow_d_ha_constants.inc, generated from the
// assembly), the reciprocal switch points (rcp14_switch.inc), and copies of the tables padded to 64
// entries for targets whose vectors hold a whole table.
#ifndef HWYPOW_CONSTANTS_H
#define HWYPOW_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define POW_CONST constexpr
#include "../svml/svml_pow_d_ha_constants.inc"
#undef POW_CONST
#include "rcp14_switch.inc"

namespace hwypow {

// rcp switch points as mantissas q in [1/2, 1): k = kmid + (q < kRcpSwitch[kmid - 32]), see hwy_pow.cc
struct RcpSwitch {
    double v[32];
    constexpr RcpSwitch() : v{} {
        constexpr unsigned p[32] = {POW_RCP14_SWITCH_PREFIXES};
        for (int i = 0; i < 32; i++) v[i] = 0.5 + p[i] * 0x1p-17;
    }
};
constexpr RcpSwitch kRcpSwitch{};

// table repeated to fill vectors of up to 64 doubles, so that TableLookupLanes can use one vector
template <size_t T>
struct alignas(64) Padded {
    double v[64];
    explicit constexpr Padded(const double *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % T];
    }
};

constexpr Padded<32> kLogThPad(vLogTh), kLogTlPad(vLogTl), kRcpSwitchPad(kRcpSwitch.v);
constexpr Padded<16> kExpThePad(vExpThe), kExpTlePad(vExpTle);

constexpr double kHalf = vHalf, kRcpLimit = vRcpLimit, kOne = vOne, kC1H = vC1H;
constexpr double kA0 = vA0, kA1 = vA1, kA2 = vA2, kA3 = vA3, kA4 = vA4, kA5 = vA5, kA6 = vA6, kA7 = vA7, kA8 = vA8,
                 kA9 = vA9;
constexpr double kB1 = vB1, kB2 = vB2, kB3 = vB3, kB4 = vB4, kB5 = vB5, kB6 = vB6;
constexpr double kShifter = vShifter, kBound = vBound;
constexpr uint64_t kRMaskBits = vRMaskBits, kExpMaskBits = vExpMaskBits;

}  // namespace hwypow

// the constant macros would clash with identifiers in Highway's headers
#undef vHalf
#undef vRcpLimit
#undef vOne
#undef vC1H
#undef vA9
#undef vA8
#undef vA7
#undef vA6
#undef vA5
#undef vA4
#undef vA3
#undef vA2
#undef vA1
#undef vA0
#undef vShifter
#undef vB6
#undef vB5
#undef vB4
#undef vB3
#undef vB2
#undef vB1
#undef vBound
#undef vRMaskBits
#undef vAbsMaskBits
#undef vExpMaskBits
#undef sL1
#undef sL2
#undef sL3
#undef sL4
#undef sE1
#undef sE2
#undef sE3
#undef sE4
#undef sE5
#undef sHuge
#undef sTiny
#undef sZero
#undef sOne
#undef sMOne
#undef sShifter
#undef sCHi
#undef sSplitter
#undef sTwo200
#undef sTwoM200
#undef POW_RCP14_SWITCH_PREFIXES

#endif
