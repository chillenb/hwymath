#!/bin/bash
# Compare the arithmetic instructions (with rounding/SAE annotations) of the compiled C vector path
# against the main path of Intel's __svml_exp8_ha, ignoring loads/moves and register allocation.
# usage: compare_insns.sh C_OBJECT ASM_OBJECT
set -euo pipefail
c_obj=$1 asm_obj=$2

ops() {  # function name, object -> sorted "mnemonic {rounding}" list of its arithmetic instructions
    objdump -d --no-show-raw-insn -M suffix "$2" \
      | awk -v f="<$1>:" '$2 == f {on = 1; next} on && (/^$/ || /\t(jne|je|ret)/) {on = 0; done = 1} on && !done' \
      | grep -oE '[[:space:]]v(fn?m(add|sub)[0-9]+pd|mulpd|subpd|addpd|scalefpd|permt2pd|permi2pd|cmppd|cmpgepd|cmpnltpd|cmp[a-z_]*pd|and(n)?pd|pandq|pandd)[[:space:]]+[^#]*' \
      | sed -E 's/^[[:space:]]+//; s/(vpermi2pd|vpermt2pd).*/vperm2pd(table)/; s/(vandpd|vpandq|vpandd).*/vand(bitwise)/; s/^vcmppd[[:space:]]+\$0x1d,/vcmpge_oqpd /; s/(vfn?m(add|sub))(132|213|231)pd/\1pd/; s/[[:space:]]+(\{[a-z-]+\})?.*/ \1/; s/ $//' \
      | sort | uniq -c
}

echo "== C (svml_exp8_ha)";      ops svml_exp8_ha "$c_obj" | tee /tmp/insns_c.$$
echo "== asm (__svml_exp8_ha)";  ops __svml_exp8_ha "$asm_obj" | tee /tmp/insns_a.$$
if diff -q /tmp/insns_c.$$ /tmp/insns_a.$$ >/dev/null; then echo "== same arithmetic instruction mix"; else echo "== DIFFERENT"; diff /tmp/insns_c.$$ /tmp/insns_a.$$ || true; fi
rm -f /tmp/insns_c.$$ /tmp/insns_a.$$
