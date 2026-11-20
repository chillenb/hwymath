"""Print C declarations for the tanh_d_ha constants, taken bit-exactly from the SVML assembly."""
import os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'common', 'tools'))
from extract_tables import extract, as_double

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_tanh_d_ha.s')
vec = extract(S, '__svml_dtanh_ha_data_internal')
sca = extract(S, '__dtanh_ha__imldTanhTab')
hx = lambda q: as_double(q).hex() if (q >> 52) & 0x7ff != 0x7ff else None

def lane(off):  # 64-byte block of 8 identical lanes
    words = vec[off // 8: off // 8 + 8]
    assert len(set(words)) == 1, off
    return words[0]

def dword(off):  # 64-byte block of 16 identical 32-bit lanes
    q = lane(off)
    assert q >> 32 == q & 0xffffffff, off
    return q & 0xffffffff

print('#ifndef TANH_TABLE  /* the Highway version includes the tables as constexpr */')
print('#define TANH_TABLE static const')
print('#endif')
print('/* vector path (__svml_dtanh_ha_data_internal): 19 tables of 16 entries, one per subinterval j of |x|')
print(' * (j = 0: |x| < 0.1875; j = 1..14: half-binades [0.1875, 0.25), [0.25, 0.375), ..., [16, 24);')
print(' * j = 15: |x| >= 24, where P0 = 1 and all other coefficients are 0).')
print(' * Row 0: B, the midpoint of the subinterval (0 for j = 0, 15), y = |x| - B; row 1: P0 (high part);')
print(' * row 2: P0 (low part);')
print(' * row k + 2: Pk, k = 1..16. tanh|x| = P0 + P1 y + ... + P16 y^16. */')
print('TANH_TABLE double vTab[19][16] = {')
for k in range(19):
    print('    {' + ', '.join(hx(vec[k * 16 + j]) for j in range(16)) + '},')
print('};')
for off, name in [(2432, 'vHiMask'), (2496, 'vHiBase'), (2560, 'vHiMax'), (10560, 'vHugeHi')]:
    print(f'#define {name:<14} 0x{dword(off):08x}')
for off, name in [(10368, 'vSignMaskBits'), (10432, 'vAbsMaskBits')]:
    print(f'#define {name:<14} 0x{lane(off):016x}ULL')

print('\n/* scalar callout (__dtanh_ha__imldTanhTab) */')
print('static const double sOnes[2] = {' + ', '.join(hx(q) for q in sca) + '};  /* +1, -1 */')
