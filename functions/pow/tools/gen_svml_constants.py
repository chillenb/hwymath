"""Print C declarations for the pow_d_ha constants, taken bit-exactly from the SVML assembly."""
import os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'common', 'tools'))
from extract_tables import extract, as_double

S = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_pow_d_ha.s')
vec = extract(S, '__svml_dpow_ha_data_internal_avx512')
sca = extract(S, '__dpow_ha_CoutTab')
assert len(vec) == 320 and len(sca) == 860
hx = lambda q: as_double(q).hex() if (q >> 52) & 0x7ff != 0x7ff else None

def lane(off):  # 64-byte block of 8 identical lanes
    words = vec[off // 8: off // 8 + 8]
    assert len(set(words)) == 1, off
    return words[0]

def table(name, words, comment, per_line=1):
    print(f'static POW_CONST double {name}[{len(words)}] = {{  /* {comment} */')
    for i in range(0, len(words), per_line):
        print('    ' + ' '.join(f'{hx(q)},' for q in words[i:i + per_line]))
    print('};')

print('#ifndef POW_CONST  /* constexpr in C++ (Highway version) */\n#define POW_CONST const\n#endif\n')
print('/* vector path (__svml_dpow_ha_data_internal_avx512) */')
table('vLogTh', vec[0:32], 'high part of -log2(rcp) (+1 for rcp >= 1.5), index = top 5 mantissa bits of rcp')
table('vLogTl', vec[32:64], 'low part')
table('vExpThe', vec[64:80], '2^(j/16), rounded')
table('vExpTle', vec[80:96], '2^(j/16) = The[j] * (1 + Tle[j])')
consts = [(768, 'vHalf'), (832, 'vRcpLimit'), (896, 'vOne'), (960, 'vC1H'),
          (1024, 'vA9'), (1088, 'vA8'), (1152, 'vA7'), (1216, 'vA6'), (1280, 'vA5'), (1344, 'vA4'),
          (1408, 'vA3'), (1472, 'vA2'), (1536, 'vA1'), (1600, 'vA0'), (1664, 'vShifter'),
          (1792, 'vB6'), (1856, 'vB5'), (1920, 'vB4'), (1984, 'vB3'), (2048, 'vB2'), (2112, 'vB1'),
          (2432, 'vBound')]
for off, name in consts:
    print(f'#define {name:<12} {hx(lane(off))}')
for off, name in [(1728, 'vRMaskBits'), (2368, 'vAbsMaskBits'), (2496, 'vExpMaskBits')]:
    print(f'#define {name:<12} 0x{lane(off):016x}ULL')

print('\n/* scalar callout (__dpow_ha_CoutTab) */')
table('sRcp1', sca[0:33], 'Rcp1 ~ 1/X1, X1 in [1,2), index = top 6 mantissa bits rounded', 3)
table('sLog1', sca[33:99], '{hi, lo} of log2(1/Rcp1)', 2)
table('sRcp2', sca[99:132], 'Rcp2 ~ 1/X2', 3)
table('sLog2', sca[132:198], '{hi, lo} of log2(1/Rcp2)', 2)
table('sRcp3', sca[198:327], 'Rcp3C ~ C/X3, C = 1/ln2', 3)
table('sLog3', sca[327:585], '{hi, lo} of log2(C/Rcp3C)', 2)
table('sExp', sca[585:841], '{hi, lo} of 2^(j/128)', 2)
for off, name in [(6728, 'sL1'), (6736, 'sL2'), (6744, 'sL3'), (6752, 'sL4'), (6760, 'sE1'), (6768, 'sE2'),
                  (6776, 'sE3'), (6784, 'sE4'), (6792, 'sE5'), (6800, 'sHuge'), (6808, 'sTiny'), (6816, 'sZero'),
                  (6824, 'sOne'), (6832, 'sMOne'), (6840, 'sShifter'), (6848, 'sCHi'), (6856, 'sSplitter'),
                  (6864, 'sTwo200'), (6872, 'sTwoM200')]:
    print(f'#define {name:<12} {hx(sca[off // 8])}')
