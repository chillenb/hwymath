// Constants for the Highway port of Intel's exp algorithm (hwy_exp_svml.cc), included once (hwy_exp_svml.cc
// itself is re-included per target). Everything comes bit-exactly from the SVML assembly
// (svml/svml_exp_d_ha_constants.inc), with the tables made constexpr so that the padded copies below are
// built at compile time.
#ifndef HWYEXPSVML_CONSTANTS_H
#define HWYEXPSVML_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define SVML_TABLE static constexpr
#include "../svml/svml_exp_d_ha_constants.inc"
#undef SVML_TABLE

namespace hwyexpsvml {

constexpr double kLog2e = vLog2e, kShifter = vShifter, kLn2Hi = vLn2Hi, kLn2Lo = vLn2Lo;
constexpr double kC1 = vC1, kC2 = vC2, kC3 = vC3, kC4 = vC4, kC5 = vC5, kC6 = vC6;
constexpr double kThreshold = vThreshold;  // |x| >= this (or NaN): scalar callout

// Table repeated to fill vectors of up to 64 doubles, so that TableLookupLanes can use one vector.
struct alignas(64) Padded16 {
    double v[64];
    explicit constexpr Padded16(const double *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % 16];
    }
};
inline const double *ThPad() { static constexpr Padded16 p(vTh); return p.v; }
inline const double *TlrPad() { static constexpr Padded16 p(vTlr); return p.v; }

}  // namespace hwyexpsvml

// the constant macros would clash with identifiers in Highway's headers (e.g. vC1)
#undef vLog2e
#undef vShifter
#undef vLn2Hi
#undef vLn2Lo
#undef vC6
#undef vC5
#undef vC4
#undef vC3
#undef vC2
#undef vC1
#undef vThreshold
#undef vRMaskBits
#undef vAbsMaskBits
#undef s64OverLn2
#undef sShifter
#undef sC2
#undef sC3
#undef sC4
#undef sC5
#undef sC6
#undef sOverflowX
#undef sDenormX
#undef sUnderflowX
#undef sLn2o64Hi
#undef sLn2o64Lo
#undef sTiny
#undef sHuge
#undef sZero
#undef sOne
#undef sTwo
#undef sTwoM60
#undef sSplitter

#endif
