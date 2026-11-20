// Constants for the Highway atan2 (hwy_atan2.cc), included once (hwy_atan2.cc itself is re-included per
// target): the SVML constant macros (bit-exact, from the assembly) turned into constexpr, because the
// macro names (vC1, ...) would clash with identifiers in Highway's headers.
#ifndef HWYATAN2_CONSTANTS_H
#define HWYATAN2_CONSTANTS_H

#include <cstdint>

#include "../svml/svml_atan2_d_ha_constants.inc"

namespace hwyatan2 {
constexpr double kAtan05Hi = vAtan05Hi, kAtan05Lo = vAtan05Lo, kPio4Hi = vPio4Hi, kPio4Lo = vPio4Lo;
constexpr double kAtan15Hi = vAtan15Hi, kAtan15Lo = vAtan15Lo, kPio2Hi = vPio2Hi, kPio2Lo = vPio2Lo;
constexpr double kHalf = vHalf, kOneHalf = vOneHalf, kOne = vOne, kPiHi = vPiHi, kPiLo = vPiLo;
constexpr double kT7 = vT7_16, kT11 = vT11_16, kT19 = vT19_16, kT39 = vT39_16;
constexpr double kC1 = vC1, kC2 = vC2, kC3 = vC3, kC4 = vC4, kC5 = vC5, kC6 = vC6, kC7 = vC7, kC8 = vC8,
                 kC9 = vC9, kC10 = vC10, kC11 = vC11, kC12 = vC12;
// the vector path handles |x|, |y| in [2^-1020, 2^993) (high word in [0x00300000, 0x7e000000))
constexpr double kMinArg = 0x1p-1020, kMaxArg = 0x1p993;
}  // namespace hwyatan2

#undef vAtan05Hi
#undef vAtan05Lo
#undef vPio4Hi
#undef vPio4Lo
#undef vAtan15Hi
#undef vAtan15Lo
#undef vPio2Hi
#undef vPio2Lo
#undef vHalf
#undef vOneHalf
#undef vPiHi
#undef vPiLo
#undef vT7_16
#undef vT19_16
#undef vT39_16
#undef vT11_16
#undef vOne
#undef vC1
#undef vC2
#undef vC3
#undef vC4
#undef vC5
#undef vC6
#undef vC7
#undef vC8
#undef vC9
#undef vC10
#undef vC11
#undef vC12
#undef vZero
#undef vZPi
#undef vZPio2
#undef vAbsMaskBits
#undef vSignMaskBits
#undef vHiBase
#undef vHiLimit

#endif
