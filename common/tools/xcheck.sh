#!/bin/bash
# Cross-architecture check of a Highway version: builds one static binary per target (x86-64 -march=native,
# AArch64 NEON/SVE/SVE2 with the conda-forge cross compiler), runs the AArch64 ones under QEMU (SVE at several
# vector lengths, which exercises every table-lookup code path: 2 lanes, 4, 8, 32), and compares the inputs and
# results bit for bit with the x86-64 run (which equals Intel's AVX-512 code: make <f>-avx512-test).
# The test program is run as `PROGRAM OUTPUT_PREFIX N` and writes OUTPUT_PREFIX.in and OUTPUT_PREFIX.out (and
# may itself check accuracy and special values; a nonzero exit status counts as a failure).
#
# usage: common/tools/xcheck.sh FUNCTION HWY_PREFIX CROSS_PREFIX [N_RANDOM]   (run from the repository root;
#        `pixi run -e cross make <f>-xcheck` does that)
# optional environment (defaults for functions/<f>/ laid out like functions/exp2):
#   XC_SRC        C++ sources        ($D/hwy/hwy_<f>.cc $D/tests/hwy_xcheck.cc)
#   XC_CSRC       C sources          ($D/svml/svml_<f>_d_ha.c: the scalar callout)
#   XC_CFLAGS     extra C flags,  XC_CXXFLAGS  extra C++ flags
#   XC_TARGETS    name:march pairs   ("neon:-march=armv8-a sve:-march=armv8.2-a+sve sve2:-march=armv9-a+sve2");
#                 names starting with "neon" run once, the others at each length in XC_VLS
#   XC_VLS        SVE vector lengths in bytes (16 32 64 256: 128, 256, 512, 2048 bits)
#   XC_SVE2_VLS   the same for names starting with "sve2" (default: XC_VLS)
#   XC_FIXED      extra fixed-length SVE builds, "name:march flags:vl" separated by ";"
set -euo pipefail
FN=$1 H=$2 E=$3 NR=${4:-2000000}
D=functions/$FN
O=build/$FN/xcheck
mkdir -p $O
SRC=${XC_SRC:-"$D/hwy/hwy_$FN.cc $D/tests/hwy_xcheck.cc"}
CSRC=${XC_CSRC:-"$D/svml/svml_${FN}_d_ha.c"}
CF="-O2 -ffp-contract=off -I. ${XC_CFLAGS:-}"
F="-O2 -ffp-contract=off -DHWY_COMPILE_ONLY_STATIC=1 -std=c++17 -I. -Icommon/include -I$H/include ${XC_CXXFLAGS:-}"
TARGETS=${XC_TARGETS:-"neon:-march=armv8-a sve:-march=armv8.2-a+sve sve2:-march=armv9-a+sve2"}
VLS=${XC_VLS:-"16 32 64 256"}
SVE2_VLS=${XC_SVE2_VLS:-$VLS}
fail=0

cobjs() {  # suffix compiler flags...: compile the C sources, print the object names
    local suf=$1 cc=$2; shift 2
    for c in $CSRC; do
        local o=$O/$(basename ${c%.c})_$suf.o
        $cc $CF "$@" -c $c -o $o
        echo $o
    done
}

# x86-64: -march=native (Highway's static AVX2 target also needs AES/PCLMUL, so x86-64-v3 is not enough)
X86_OBJS=$(cobjs x86 gcc -march=native)
g++ $F -march=native $SRC $X86_OBJS -lm -o $O/xcheck_x86
./$O/xcheck_x86 $O/x86 "$NR" || fail=1

CC=$E/bin/aarch64-conda-linux-gnu-gcc CXX=$E/bin/aarch64-conda-linux-gnu-g++ QEMU=$E/bin/qemu-aarch64
A64_OBJS=$(cobjs a64 $CC)
check() {  # name
    cmp -s $O/x86.in $O/$1.in && cmp -s $O/x86.out $O/$1.out && echo "    -> $1: bit-identical to x86-64" \
        || { echo "    -> $1: DIFFERENT"; fail=1; }
}
run() {  # name qemu-cpu-options
    $QEMU -cpu $2 ./$O/xcheck_a64_${1%%@*} $O/a64_${1/@/_} "$NR" || fail=1
    check a64_${1/@/_}
}
for v in $TARGETS; do
    name=${v%%:*} march=${v#*:}
    $CXX $F $march -static $SRC $A64_OBJS -lm -o $O/xcheck_a64_$name
    case $name in
        neon*) run $name max ;;
        sve2*) for vl in $SVE2_VLS; do run $name@$vl max,sve-default-vector-length=$vl; done ;;
        *)     for vl in $VLS; do run $name@$vl max,sve-default-vector-length=$vl; done ;;
    esac
done
if [ -n "${XC_FIXED:-}" ]; then
    IFS=';' read -ra fixed <<< "$XC_FIXED"
    for v in "${fixed[@]}"; do
        IFS=: read -r name march vl <<< "$v"
        $CXX $F $march -static $SRC $A64_OBJS -lm -o $O/xcheck_a64_$name
        run $name max,sve-default-vector-length=$vl
    done
fi
[ $fail = 0 ] && echo "$FN xcheck: all checks passed" || echo "$FN xcheck: FAILED"
exit $fail
