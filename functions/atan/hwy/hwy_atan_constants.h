// Constants for the Highway atan (hwy_atan.cc), included once (hwy_atan.cc itself is re-included per
// target): the SVML constants (bit-exact, from the assembly) as constexpr, because the macro names would
// clash with identifiers in Highway's headers, and the 32-entry tables repeated to fill vectors of up to
// 64 doubles (compile-time constants, so that no run-time initialization ends up in the loop).
#ifndef HWYATAN_CONSTANTS_H
#define HWYATAN_CONSTANTS_H

#include <cstddef>
#include <cstdint>

#define static static constexpr  // the tables in the .inc become constexpr
#include "../svml/svml_atan_d_ha_constants.inc"
#undef static

namespace hwyatan {
constexpr double kShifter = vShifter, kTabMax = vTabMax, kMOne = vMOne, kOne = vOne, kYMax = vYMax;
constexpr double kPio2Hi = vPio2Hi, kPio2Lo = vPio2Lo;
constexpr double kA0 = vA0, kA1 = vA1, kA2 = vA2, kA3 = vA3, kA4 = vA4, kA5 = vA5;

struct alignas(64) Padded {  // 32-entry table repeated to 64 entries (TableLookupLanes on one vector)
    double v[64];
    explicit constexpr Padded(const double *t) : v{} {
        for (size_t i = 0; i < 64; i++) v[i] = t[i % 32];
    }
};
inline const double *AtanHi() { return vAtanHi; }
inline const double *AtanLo() { return vAtanLo; }
inline const double *AtanHiPad() { static constexpr Padded p(vAtanHi); return p.v; }
inline const double *AtanLoPad() { static constexpr Padded p(vAtanLo); return p.v; }
}  // namespace hwyatan

#undef vShifter
#undef vTabMax
#undef vMOne
#undef vOne
#undef vYMax
#undef vShifter16
#undef vPio2Hi
#undef vPio2Lo
#undef vA5
#undef vA4
#undef vA3
#undef vA2
#undef vA1
#undef vA0
#undef vAbsMaskBits

#endif
