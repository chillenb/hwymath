"""Cross-check of the Highway exp results collected by CI (`hwy_test --hash FILE` on several machines).

usage: python3 tools/check_hashes.py HASH_FILE...

Each line is "arch target variant n hash" (FNV-1a of the result bits) or "arch inputs - n hash".
Checks that all runs used the same inputs, and that every target running the same hp variant
(hp4: 4-entry table, hp8: 8-entry table) produced bit-identical results on every machine.
Exit status 1 on any mismatch, or if no hp4 results from both x86_64 and aarch64 were found.
"""
import collections
import sys


def main(paths):
    inputs, results = set(), collections.defaultdict(set)
    rows = []
    for path in paths:
        for line in open(path):
            arch, target, variant, n, h = line.split()
            if target == 'inputs':
                inputs.add((n, h))
            else:
                results[(variant, n)].add(h)
                rows.append((variant, arch, target, h, path))
    ok = True
    print(f'inputs: {len(inputs)} distinct input set(s) {sorted(inputs)}')
    if len(inputs) != 1:
        print('FAIL: the runs used different inputs')
        ok = False
    for variant, arch, target, h, path in sorted(rows):
        print(f'  {variant}  {arch:8} {target:14} {h}  ({path})')
    for (variant, n), hashes in sorted(results.items()):
        status = 'identical' if len(hashes) == 1 else 'DIFFERENT'
        print(f'{variant} (n={n}): {status} across all machines and targets')
        ok &= len(hashes) == 1
    archs_hp4 = {arch for variant, arch, *_ in rows if variant == 'hp4'}
    if not {'x86_64', 'aarch64'} <= archs_hp4:
        print(f'FAIL: expected hp4 results from x86_64 and aarch64, got {sorted(archs_hp4)}')
        ok = False
    print('cross-check passed' if ok else 'cross-check FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
