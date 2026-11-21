#!/bin/bash
# Cross-architecture check of the two Highway exps (hp, and svml = Intel's algorithm): builds one static
# binary per target (x86-64 AVX2 natively, AArch64 NEON/SVE/SVE2 with the conda-forge cross compiler),
# runs the AArch64 ones under QEMU (SVE at 128/256/512/2048-bit vectors: 2, 4, 8 and 32 lanes, i.e.
# every table-lookup variant), and compares inputs and both results bit for bit with the x86-64 run
# (whose svml results equal the C translation: exp-test).
# usage: functions/exp/tools/xcheck.sh HWY_PREFIX CROSS_PREFIX [N_RANDOM]   (run from the repository root;
#        `pixi run -e cross make exp-xcheck` does that)
set -euo pipefail
H=$1 E=$2 NR=${3:-2000000}
E_DIR=functions/exp
O=build/exp/xcheck
mkdir -p $O
F="-O2 -ffp-contract=off -DHWY_COMPILE_ONLY_STATIC=1 -std=c++17 -I. -Icommon/include -I$H/include"
SRC="$E_DIR/hwy/hwy_exp.cc $E_DIR/hwy/hwy_exp_svml.cc $E_DIR/tests/hwy_xcheck.cc"

# x86-64: -march=native (Highway's static AVX2 target also needs AES/PCLMUL, so x86-64-v3 is not enough)
gcc -O2 -ffp-contract=off -march=native -c $E_DIR/svml/svml_exp_d_ha.c -o $O/rare_x86.o
g++ $F -march=native $SRC $O/rare_x86.o -lm -o $O/xcheck_x86
./$O/xcheck_x86 $O/x86 "$NR"

CC=$E/bin/aarch64-conda-linux-gnu-gcc CXX=$E/bin/aarch64-conda-linux-gnu-g++ QEMU=$E/bin/qemu-aarch64
$CC -O2 -ffp-contract=off -c $E_DIR/svml/svml_exp_d_ha.c -o $O/rare_a64.o
fail=0
check() {  # name
    cmp -s $O/x86.in $O/$1.in || { echo "    -> $1: DIFFERENT inputs"; fail=1; return; }
    for r in out:hp svml.out:svml; do
        cmp -s $O/x86.${r%%:*} $O/$1.${r%%:*} && echo "    -> $1 ${r#*:}: bit-identical to x86-64" \
            || { echo "    -> $1 ${r#*:}: DIFFERENT"; fail=1; }
    done
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
