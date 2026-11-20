#!/bin/bash
# Cross-architecture check of the Highway pow: builds one static binary per target (x86-64 AVX2 natively,
# AArch64 NEON/SVE/SVE2 with the conda-forge cross compiler), runs the AArch64 ones under QEMU (SVE at
# 128/256/512/2048-bit vectors), and compares inputs and results bit for bit with the x86-64 run (whose AVX2
# results equal Intel's AVX-512 code: make pow-avx512-test).
# usage: functions/pow/tools/xcheck.sh HWY_PREFIX CROSS_PREFIX [N_RANDOM]   (run from the repository root;
#        `pixi run -e cross make pow-xcheck` does that)
set -euo pipefail
H=$1 E=$2 NR=${3:-1000000}
P=functions/pow
O=build/pow/xcheck
mkdir -p $O
F="-O2 -ffp-contract=off -DHWY_COMPILE_ONLY_STATIC=1 -std=c++17 -I. -Icommon/include -I$H/include"
SRC="$P/hwy/hwy_pow.cc $P/tests/hwy_xcheck.cc"

# x86-64: -march=native (Highway's static AVX2 target also needs AES/PCLMUL, so x86-64-v3 is not enough)
gcc -O2 -ffp-contract=off -march=native -c $P/svml/svml_pow_d_ha.c -o $O/rare_x86.o
g++ $F -march=native $SRC $O/rare_x86.o -lm -o $O/xcheck_x86
./$O/xcheck_x86 $O/x86 "$NR"

CC=$E/bin/aarch64-conda-linux-gnu-gcc CXX=$E/bin/aarch64-conda-linux-gnu-g++ QEMU=$E/bin/qemu-aarch64
$CC -O2 -ffp-contract=off -c $P/svml/svml_pow_d_ha.c -o $O/rare_a64.o
fail=0
check() {  # name
    cmp -s $O/x86.in $O/$1.in && cmp -s $O/x86.out $O/$1.out && echo "    -> $1: bit-identical to x86-64" || { echo "    -> $1: DIFFERENT"; fail=1; }
}
for v in "neon:-march=armv8-a" "sve:-march=armv8.2-a+sve" "sve2:-march=armv9-a+sve2"; do
    name=${v%%:*} march=${v#*:}
    $CXX $F $march -static $SRC $O/rare_a64.o -lm -o $O/xcheck_a64_$name
    if [ $name = neon ]; then
        $QEMU -cpu max ./$O/xcheck_a64_$name $O/a64_$name "$NR"; check a64_$name
    else
        for vl in 16 32 64 256; do  # vector length in bytes: 128, 256, 512, 2048 bits
            $QEMU -cpu max,sve-default-vector-length=$vl ./$O/xcheck_a64_$name $O/a64_${name}_$vl "$NR"
            check a64_${name}_$vl
        done
    fi
done
exit $fail
