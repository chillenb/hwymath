// Constants of the Highway log1p (hwy_log1p.cc), included once (hwy_log1p.cc itself is re-included per target).
// hwy_log1p_constants.inc is generated from the SVML assembly by functions/log/tools/gen_svml_constants.py log1p --cxx.
#ifndef HWYMATH_HWY_LOG1P_CONSTANTS_H
#define HWYMATH_HWY_LOG1P_CONSTANTS_H

#include <cstdint>

#include "functions/log/hwy/log_tables.h"

namespace hwymath_log {
namespace log1p_e {
#include "hwy_log1p_constants.inc"
static const double *ThPad() { static constexpr Padded16 p(kTh); return p.v; }
static const double *TlPad() { static constexpr Padded16 p(kTl); return p.v; }
}  // namespace log1p_e
}  // namespace hwymath_log

#endif
