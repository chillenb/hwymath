// Constants of the Highway log (hwy_log.cc), included once (hwy_log.cc itself is re-included per target).
// hwy_log_constants.inc is generated from the SVML assembly by tools/gen_svml_constants.py log --cxx.
#ifndef HWYMATH_HWY_LOG_CONSTANTS_H
#define HWYMATH_HWY_LOG_CONSTANTS_H

#include <cstdint>

#include "functions/log/hwy/log_tables.h"

namespace hwymath_log {
namespace log_e {
#include "hwy_log_constants.inc"
static const double *ThPad() { static constexpr Padded16 p(kTh); return p.v; }
static const double *TlPad() { static constexpr Padded16 p(kTl); return p.v; }
}  // namespace log_e
}  // namespace hwymath_log

#endif
