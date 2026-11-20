"""Print C declarations for the cosh_d_ha constants, taken bit-exactly from the SVML assembly."""
import os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'common', 'tools'))
from extract_tables import extract, as_double

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_cosh_d_ha.s')
vec = extract(S, '__svml_dcosh_ha_data_internal')
sca = extract(S, '__dcosh_ha_CoutTab')
hx = lambda q: as_double(q).hex() if (q >> 52) & 0x7ff != 0x7ff else None

def lane(off):  # 64-byte block of 8 identical lanes
    words = vec[off // 8: off // 8 + 8]
    assert len(set(words)) == 1, off
    return words[0]

def table(name, off, comment):
    print(f'COSH_TABLE double {name}[16] = {{  /* {comment} */')
    for j in range(16): print(f'    {hx(vec[off // 8 + j])},')
    print('};')

print('#ifndef COSH_TABLE  /* the Highway version includes the tables as constexpr */')
print('#define COSH_TABLE static const')
print('#endif')
print('/* vector path (__svml_dcosh_ha_data_internal); j = 0..15 */')
table('vTpH', 0, '2^(j/16)/2, rounded')
table('vTpL', 128, '2^(j/16)/2 - vTpH[j]')
table('vTmH', 256, '2^(-j/16)/2, rounded')
for off, name in [(384, 'vShifter'), (512, 'vC2'), (576, 'vS3'), (640, 'vC4'), (704, 'vS5'), (768, 'vC6'),
                  (832, 'vS7'), (3968, 'vLog2e'), (4032, 'vLn2Hi'), (4096, 'vLn2Lo')]:
    print(f'#define {name:<14} {hx(lane(off))}')
for off, name in [(4608, 'vExpMaskBits'), (4672, 'vSignMaskBits')]:
    print(f'#define {name:<14} 0x{lane(off):016x}ULL')
thr = lane(4736)
assert thr >> 32 == thr & 0xffffffff
print(f'#define {"vThresholdHi":<14} 0x{thr & 0xffffffff:08x}  /* callout if the high word of |x| is greater (signed) */')

print('\n/* scalar callout (__dcosh_ha_CoutTab) */')
print('static const double sTab[130] = {  /* {hi, lo} of 2^(j/64), j = 0..63; then 2.0, 0.0 (j = 64) */')
for j in range(65): print(f'    {hx(sca[2*j])}, {hx(sca[2*j+1])},')
print('};')
for off, name in [(1024, 'sTwo'), (1040, 's64OverLn2'), (1048, 'sShifter'), (1056, 'sP2'), (1064, 'sP3'), (1072, 'sP4'),
                  (1080, 'sP5'), (1088, 'sP6'), (1096, 'sOverflowX'), (1104, 'sLn2o64Hi'), (1112, 'sLn2o64Lo'),
                  (1120, 'sHuge'), (1136, 'sOne'), (1144, 'sLargeX')]:
    print(f'#define {name:<14} {hx(sca[off // 8])}')
assert sca[1128 // 8] == 0
