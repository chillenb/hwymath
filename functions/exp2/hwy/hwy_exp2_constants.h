// Constants for the Highway exp2 (hwy_exp2.cc), included once (hwy_exp2.cc itself is re-included per
// target). Everything comes bit-exactly from the SVML assembly (svml/svml_exp2_d_ha_constants.inc), with
// the tables made constexpr so that the padded copies below are built at compile time.
#ifndef HWYEXP2_CONSTANTS_H
#define HWYEXP2_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define SVML_TABLE static constexpr
#include "../svml/svml_exp2_d_ha_constants.inc"
#undef SVML_TABLE

namespace hwyexp2 {

constexpr double kC1 = vC1, kC2 = vC2, kC3 = vC3, kC4 = vC4, kC5 = vC5, kC6 = vC6;
constexpr double kThreshold = vThreshold;  // |x| >= this (or inf): scalar callout
// 1.5*2^52 + 1023*16: Floor(16 x) + this holds floor(16 x) + 16*1023 in its low mantissa bits, so its
// low 4 bits are the table index and bits 4..15 the biased exponent of 2^floor(x)
constexpr double kShifter52 = 0x1.8p52 + 1023 * 16;
constexpr double kShifter = vShifter;  // vaddpd {rd-sae} x + kShifter: floor(16 x) in the low bits (AVX-512 path)

// Table repeated to fill vectors of up to 64 doubles, so that TableLookupLanes can use one vector.
struct alignas(64) Padded16 {
    double v[64];
    explicit constexpr Padded16(const double *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % 16];
    }
};
inline const double *ThPad() { static constexpr Padded16 p(vTh); return p.v; }
inline const double *TlPad() { static constexpr Padded16 p(vTl); return p.v; }

}  // namespace hwyexp2

// the constant macros would clash with identifiers in Highway's headers (e.g. vC1)
#undef vC1
#undef vC2
#undef vC3
#undef vC4
#undef vC5
#undef vC6
#undef vShifter
#undef vThreshold
#undef vAbsMaskBits
#undef vIndexMaskBits
#undef sShifter
#undef sC1
#undef sC2
#undef sC3
#undef sC4
#undef sC5
#undef sOverflowX
#undef sDenormX
#undef sUnderflowX
#undef sTiny
#undef sHuge
#undef sZero
#undef sTwo
#undef sTwoM100

#endif
