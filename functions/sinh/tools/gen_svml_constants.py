"""Print C declarations for the sinh_d_ha constants, taken bit-exactly from the SVML assembly."""
import os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'common', 'tools'))
from extract_tables import extract, as_double

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_sinh_d_ha.s')
vec = extract(S, '__svml_dsinh_ha_data_internal')
sca = extract(S, '__dsinh_ha_CoutTab')
hx = lambda q: as_double(q).hex() if (q >> 52) & 0x7ff != 0x7ff else None

def lane(off):  # 64-byte block of 8 identical lanes
    words = vec[off // 8: off // 8 + 8]
    assert len(set(words)) == 1, off
    return words[0]

def table(name, off, comment):
    print(f'SINH_TABLE double {name}[16] = {{  /* {comment} */')
    for j in range(16): print(f'    {hx(vec[off // 8 + j])},')
    print('};')

print('#ifndef SINH_TABLE  /* the Highway version includes the tables as constexpr */')
print('#define SINH_TABLE static const')
print('#endif')
print('/* vector path (__svml_dsinh_ha_data_internal); j = 0..15 */')
table('vTpH', 0, '2^(j/16)/2, rounded')
table('vTpL', 128, '2^(j/16)/2 - vTpH[j]')
table('vTmH', 256, '2^(-j/16)/2, rounded')
table('vTmL', 384, '2^(-j/16)/2 - vTmH[j]')
for off, name in [(512, 'vShifter'), (768, 'vC2'), (832, 'vS3'), (896, 'vC4'), (960, 'vS5'), (1024, 'vC6'),
                  (1088, 'vS7'), (1152, 'vC8'), (1216, 'vLog2e'), (1280, 'vLn2Hi'), (1344, 'vLn2Lo')]:
    print(f'#define {name:<14} {hx(lane(off))}')
for off, name in [(576, 'vShifterP1Bits'), (704, 'vIndexMaskBits'), (1408, 'vSignMaskBits')]:
    print(f'#define {name:<14} 0x{lane(off):016x}ULL')
thr = lane(640)
assert thr >> 32 == thr & 0xffffffff
print(f'#define {"vThresholdHi":<14} 0x{thr & 0xffffffff:08x}  /* callout if the high word of |x| is greater (signed) */')

print('\n/* scalar callout (__dsinh_ha_CoutTab) */')
print('static const double sTab[130] = {  /* {hi, lo} of 2^(j/64), j = 0..63; then 2.0, 0.0 (j = 64) */')
for j in range(65): print(f'    {hx(sca[2*j])}, {hx(sca[2*j+1])},')
print('};')
for off, name in [(1024, 'sTwo'), (1040, 'sP2'), (1048, 'sP3'), (1056, 'sP4'), (1064, 'sP5'), (1072, 'sP6'),
                  (1080, 'sS3'), (1088, 'sS5'), (1096, 'sS7'), (1104, 'sS9'), (1112, 's64OverLn2'), (1120, 'sShifter'),
                  (1128, 'sSplitter'), (1136, 'sOverflowX'), (1144, 'sLn2o64Hi'), (1152, 'sLn2o64Lo'),
                  (1160, 'sTiny'), (1168, 'sHuge'), (1176, 'sSmallX'), (1184, 'sLargeX')]:
    print(f'#define {name:<14} {hx(sca[off // 8])}')
