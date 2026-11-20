#!/bin/bash
# Compare the arithmetic instructions (mnemonic + rounding/SAE annotation) of a compiled C vector function
# with the main path of an SVML assembly routine, ignoring loads/moves, register allocation and the
# choice of FMA form (132/213/231). Works for double, single and half precision (pd/ps/ph).
# The main path ends at the first conditional jump or return (the branch to the scalar callout).
#
# usage: compare_insns.sh C_OBJECT C_SYMBOL ASM_OBJECT ASM_SYMBOL
#   e.g. compare_insns.sh build/exp/svml_exp_d_ha.o svml_exp8_ha build/exp/svml_asm.o __svml_exp8_ha
set -euo pipefail
[ $# -eq 4 ] || { sed -n '2,8p' "$0"; exit 2; }
c_obj=$1 c_sym=$2 asm_obj=$3 asm_sym=$4

ops() {  # symbol, object -> sorted "mnemonic {rounding}" counts of its arithmetic instructions
    objdump -d --no-show-raw-insn -M suffix "$2" \
      | awk -v f="<$1>:" '$2 == f {on = 1; next} on && (/^$/ || /\t(jne|je|ret)/) {on = 0; done = 1} on && !done' \
      | grep -oE '[[:space:]]v(fn?m(add|sub)[0-9]+p[sdh]|(mul|sub|add|div|min|max|sqrt|scalef|rcp14|rsqrt14|getexp|getmant|rndscale|reduce|fixupimm|range)p[sdh]|perm[it]2p[sdh]|permp[sdh]|cmp[a-z_]*p[sdh]|(and|andn|or|xor)p[sdh]|p(and|andn|or|xor)[dq])[[:space:]]+[^#]*' \
      | sed -E 's/^[[:space:]]+//
                s/(vperm[it]2p[sdh]).*/vperm2(table)/
                s/(v(and|andn|or|xor)p[sdh]|vp(and|andn|or|xor)[dq])[[:space:]].*/v\2\3(bitwise)/
                s/^vcmpp([sdh])[[:space:]]+\$0x1d,/vcmpge_oqp\1 /
                s/(vfn?m(add|sub))(132|213|231)p([sdh])/\1p\4/
                s/[[:space:]]+(\{[a-z-]+\})?.*/ \1/
                s/ $//' \
      | sort | uniq -c
}

tc=$(mktemp) ta=$(mktemp)
trap 'rm -f "$tc" "$ta"' EXIT
echo "== C ($c_sym)";     ops "$c_sym" "$c_obj" | tee "$tc"
echo "== asm ($asm_sym)"; ops "$asm_sym" "$asm_obj" | tee "$ta"
if diff -q "$tc" "$ta" >/dev/null; then echo "== same arithmetic instruction mix"; else echo "== DIFFERENT"; diff "$tc" "$ta" || true; exit 1; fi
