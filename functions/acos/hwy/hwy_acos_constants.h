// Constants for the Highway acos (hwy_acos.cc), included once (hwy_acos.cc itself is re-included per
// target): the SVML constant macros (bit-exact, from the assembly) turned into constexpr, because the
// macro names (vC1, ...) would clash with identifiers in Highway's headers.
#ifndef HWYACOS_CONSTANTS_H
#define HWYACOS_CONSTANTS_H

#include <cstdint>

#include "../svml/svml_acos_d_ha_constants.inc"

namespace hwyacos {
constexpr double kHalf = vHalf, kSmallY = vSmallY, kMOne = vMOne, kTwo = vTwo;
constexpr double kSqrtC0 = vSqrtC0, kSqrtC1 = vSqrtC1, kSqrtC2 = vSqrtC2, kSqrtC3 = vSqrtC3;
constexpr double kC1 = vC1, kC2 = vC2, kC3 = vC3, kC4 = vC4, kC5 = vC5, kC6 = vC6, kC7 = vC7, kC8 = vC8,
                 kC9 = vC9, kC10 = vC10, kC11 = vC11, kC12 = vC12;
constexpr double kPio2Hi = vPio2Hi, kPio2Lo = vPio2Lo, kPiHi = vPiHi, kPiLo = vPiLo;
}  // namespace hwyacos

#undef vHalf
#undef vSmallY
#undef vMOne
#undef vTwo
#undef vSqrtC3
#undef vSqrtC2
#undef vSqrtC1
#undef vSqrtC0
#undef vC12
#undef vC11
#undef vC10
#undef vC9
#undef vC8
#undef vC7
#undef vC6
#undef vC5
#undef vC4
#undef vC3
#undef vC2
#undef vC1
#undef vPio2Hi
#undef vPio2Lo
#undef vZero
#undef vSignMaskBits
#undef vPiLoBits
#undef vPiHiBits
#undef vPiLo
#undef vPiHi

#endif
