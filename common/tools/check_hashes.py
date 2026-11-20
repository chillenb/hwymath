"""Cross-check of result hashes collected by CI (`make hashes` on several machines).

usage: python3 common/tools/check_hashes.py HASH_FILE...

Each line is "arch function target variant n hash" (FNV-1a of the result bits), or
"arch function inputs - n hash" for the inputs. For every function, checks that all runs used the same
inputs, that all targets with the same variant (e.g. exp's hp4/hp8: hp algorithm with a 4/8-entry
table) produced bit-identical results on every machine, and that at least one variant was compared
across x86_64 and aarch64. Exit status 1 on any failure.
"""
import collections
import sys


def main(paths):
    inputs = collections.defaultdict(set)    # function -> {(n, hash)}
    results = collections.defaultdict(set)   # (function, variant, n) -> {hash}
    archs = collections.defaultdict(set)     # (function, variant) -> {arch}
    rows = []
    for path in paths:
        for line in open(path):
            if not line.strip():
                continue
            arch, func, target, variant, n, h = line.split()
            if target == 'inputs':
                inputs[func].add((n, h))
            else:
                results[(func, variant, n)].add(h)
                archs[(func, variant)].add(arch)
                rows.append((func, variant, arch, target, h, path))
    ok = True
    for func in sorted(set(inputs) | {f for f, *_ in results}):
        print(f'== {func}')
        if len(inputs[func]) != 1:
            print(f'  FAIL: {len(inputs[func])} different input sets {sorted(inputs[func])}')
            ok = False
        for f, variant, arch, target, h, path in sorted(rows):
            if f == func:
                print(f'  {variant:5} {arch:8} {target:18} {h}  ({path})')
        for (f, variant, n), hashes in sorted(results.items()):
            if f == func:
                same = len(hashes) == 1
                print(f'  {variant} (n={n}): {"identical" if same else "DIFFERENT"} on {sorted(archs[(f, variant)])}')
                ok &= same
        if not any({'x86_64', 'aarch64'} <= a for (f, _), a in archs.items() if f == func):
            print('  FAIL: no variant was compared between x86_64 and aarch64')
            ok = False
    print('cross-check passed' if ok else 'cross-check FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
