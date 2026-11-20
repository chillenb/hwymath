// Constants of the Highway log10 (hwy_log10.cc), included once (hwy_log10.cc itself is re-included per target).
// hwy_log10_constants.inc is generated from the SVML assembly by functions/log/tools/gen_svml_constants.py log10 --cxx.
#ifndef HWYMATH_HWY_LOG10_CONSTANTS_H
#define HWYMATH_HWY_LOG10_CONSTANTS_H

#include <cstdint>

#include "functions/log/hwy/log_tables.h"

namespace hwymath_log {
namespace log10_e {
#include "hwy_log10_constants.inc"
static const double *ThPad() { static constexpr Padded16 p(kTh); return p.v; }
static const double *TlPad() { static constexpr Padded16 p(kTl); return p.v; }
}  // namespace log10_e
}  // namespace hwymath_log

#endif
