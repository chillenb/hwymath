"""Print C declarations for the exp_d_ha constants, taken bit-exactly from the SVML assembly."""
import os, struct, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_tables import extract, as_double

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_exp_d_ha.s')
vec = extract(S, '__svml_dexp_ha_data_internal_avx512')
sca = extract(S, '_imldExpHATab')
hx = lambda q: as_double(q).hex() if (q >> 52) & 0x7ff != 0x7ff else None

def lane(off):  # 64-byte block of 8 identical lanes
    words = vec[off // 8: off // 8 + 8]
    assert len(set(words)) == 1, off
    return words[0]

print('/* vector path (__svml_dexp_ha_data_internal_avx512) */')
print('static const double vTh[16] = {  /* 2^(j/16), rounded */')
for j in range(16): print(f'    {hx(vec[j])},')
print('};')
print('static const double vTlr[16] = {  /* 2^(j/16) = Th[j] * (1 + Tlr[j]) */')
for j in range(16): print(f'    {hx(vec[16 + j])},')
print('};')
for off, name in [(256, 'vLog2e'), (320, 'vShifter'), (384, 'vLn2Hi'), (448, 'vLn2Lo'), (576, 'vC6'), (640, 'vC5'),
                  (704, 'vC4'), (768, 'vC3'), (832, 'vC2'), (896, 'vC1'), (1024, 'vThreshold')]:
    print(f'#define {name:<12} {hx(lane(off))}')
for off, name in [(512, 'vRMaskBits'), (960, 'vAbsMaskBits')]:
    print(f'#define {name:<12} 0x{lane(off):016x}ULL')

print('\n/* scalar callout (_imldExpHATab) */')
print('static const double sTab[128] = {  /* {hi, lo} of 2^(j/64), j = 0..63 */')
for j in range(64): print(f'    {hx(sca[2*j])}, {hx(sca[2*j+1])},')
print('};')
for off, name in [(1024, 's64OverLn2'), (1032, 'sShifter'), (1040, 'sC2'), (1048, 'sC3'), (1056, 'sC4'), (1064, 'sC5'),
                  (1072, 'sC6'), (1080, 'sOverflowX'), (1088, 'sDenormX'), (1096, 'sUnderflowX'), (1104, 'sLn2o64Hi'),
                  (1112, 'sLn2o64Lo'), (1120, 'sTiny'), (1128, 'sHuge'), (1136, 'sZero'), (1144, 'sOne'),
                  (1152, 'sTwo'), (1160, 'sTwoM60'), (1168, 'sSplitter')]:
    print(f'#define {name:<12} {hx(sca[off // 8])}')
