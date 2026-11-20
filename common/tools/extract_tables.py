"""Extract a data object from an SVML .s file as a list of 64-bit words (exact bit patterns)."""
import re, struct, sys

def extract(path, symbol):
    src = open(path).read()
    start = src.index(f'\n{symbol}:') + len(symbol) + 2
    words32 = []
    for line in src[start:].splitlines():
        m = re.match(r'\s*\.(long|quad)\s+(.*)', line)
        if m:
            for v in m.group(2).split(','):
                v = int(v.strip(), 0)
                words32 += [v] if m.group(1) == 'long' else [v & 0xffffffff, v >> 32]
        elif re.match(r'\s*\.(align|byte)', line) or not line.strip():
            continue
        else:
            break
    raw = struct.pack('<%dI' % len(words32), *words32)
    return [struct.unpack_from('<Q', raw, i)[0] for i in range(0, len(raw), 8)]

def as_double(q):
    return struct.unpack('<d', struct.pack('<Q', q))[0]

if __name__ == '__main__':
    for q in extract(sys.argv[1], sys.argv[2]):
        print(f'0x{q:016x}  {as_double(q)!r}')
