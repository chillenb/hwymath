"""Mechanical translation of an SVML scalar callout routine (x86-64 SSE2 scalar code) into portable C.

usage: asm2c.py ASM_FILE ASM_SYMBOL C_NAME NARGS > out.inc

Each instruction becomes one C statement on a small machine model (shown next to it as a comment), so
the C function computes exactly what the assembly computes, operation by operation:
  - xmm registers hold doubles (only the low lane is ever used by these routines); scalar SSE2
    arithmetic is IEEE binary64 with round-to-nearest, i.e. plain C double arithmetic (compile without
    FMA contraction);
  - general-purpose registers are uint64_t with x86 sub-register semantics (32-bit writes zero-extend,
    8/16-bit writes merge);
  - memory is real memory: a local stack array for the red zone below %rsp, the caller's pointers in
    %rdi/%rsi/%rdx, and the routine's data tables (emitted as uint64_t arrays with the assembly's bits);
  - condition flags are kept as the operands of the last flag-setting instruction and evaluated by
    the conditional jumps; comisd uses signaling comparisons and ucomisd quiet ones (isgreater, ...),
    so the same floating-point exceptions are raised; little-endian byte order is assumed.
The tables are writable (not const) so that the compiler cannot constant-fold operations on table
values (which would drop the floating-point exceptions the assembly raises at run time); zeroing a
register (pxor) reads a volatile zero for the same reason. GCC/Clang optimize the model into ordinary
register code.
"""
import os, re, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')
sys.path.insert(0, os.path.join(ROOT, 'common', 'tools'))
from extract_tables import extract  # noqa: E402

asm_path, sym, cname, nargs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
src = open(asm_path).read().split('\n')

# --- the routine's lines ---
start = src.index(f'{sym}:')
body = []
for line in src[start + 1:]:
    if re.match(r'\s*\.type\s', line) or re.match(r'\s*\.size\s', line):
        break
    body.append(line)

GPR64 = ['rax', 'rbx', 'rcx', 'rdx', 'rsi', 'rdi', 'rbp', 'rsp'] + [f'r{i}' for i in range(8, 16)]
ALIAS = {}
for r in ['ax', 'bx', 'cx', 'dx']:
    ALIAS['r' + r] = ('r' + r, 64); ALIAS['e' + r] = ('r' + r, 32); ALIAS[r] = ('r' + r, 16)
    ALIAS[r[0] + 'l'] = ('r' + r, 8)
for r in ['si', 'di', 'bp', 'sp']:
    ALIAS['r' + r] = ('r' + r, 64); ALIAS['e' + r] = ('r' + r, 32); ALIAS[r] = ('r' + r, 16)
    ALIAS[r + 'l'] = ('r' + r, 8)
for i in range(8, 16):
    ALIAS[f'r{i}'] = (f'r{i}', 64); ALIAS[f'r{i}d'] = (f'r{i}', 32); ALIAS[f'r{i}w'] = (f'r{i}', 16)
    ALIAS[f'r{i}b'] = (f'r{i}', 8)
SFX = {'b': 8, 'w': 16, 'l': 32, 'q': 64}
UT = {8: 'uint8_t', 16: 'uint16_t', 32: 'uint32_t', 64: 'uint64_t'}
ST = {8: 'int8_t', 16: 'int16_t', 32: 'int32_t', 64: 'int64_t'}

tables = []  # data symbols referenced


def tab_name(s):
    if s not in tables:
        tables.append(s)
    return 'T_' + re.sub(r'[^A-Za-z0-9]', '_', s).strip('_')


def split_ops(s):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur.strip()); cur = ''
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def addr(op):
    """C expression (uint64_t) for a memory operand."""
    m = re.fullmatch(r'(-?\d+)?\+?([._A-Za-z][._A-Za-z0-9]*)?\(%rip\)', op)
    if m:
        off = int(m.group(1) or 0)
        return f'((uint64_t)(uintptr_t){tab_name(m.group(2))} + {off})'
    m = re.fullmatch(r'(-?\d+)?\(%(\w+)(?:,%(\w+),(\d+))?\)', op)
    assert m, op
    off, base, idx, sc = int(m.group(1) or 0), m.group(2), m.group(3), m.group(4)
    e = f'{ALIAS[base][0]}'
    if idx:
        e += f' + {ALIAS[idx][0]} * {sc}'
    return f'({e} + {off})' if off else f'({e})'


def is_mem(op):
    return '(' in op


def is_xmm(op):
    return op.startswith('%xmm')


def x(op):
    return 'x' + op[4:]


def rd(op, w):
    """read an integer operand of width w"""
    if op.startswith('$'):
        return f'(({UT[w]})({int(op[1:], 0)}))'
    if is_mem(op):
        return f'ld{w}({addr(op)})'
    r, rw = ALIAS[op[1:]]
    assert rw == w, (op, w)
    return f'(({UT[w]}){r})'


def wr(op, w, val):
    """write an integer operand of width w"""
    if is_mem(op):
        return f'st{w}({addr(op)}, {val});'
    r, rw = ALIAS[op[1:]]
    assert rw == w, (op, w)
    if w == 64:
        return f'{r} = (uint64_t)({val});'
    if w == 32:
        return f'{r} = (uint32_t)({val});'
    mask = (1 << w) - 1
    return f'{r} = ({r} & ~(uint64_t){mask:#x}) | ({UT[w]})({val});'


def flags_int(res, w, a='', b=''):
    if a:
        return f'fk = FK_INT; fw = {w}; fa = (uint64_t)({a}); fb = (uint64_t)({b});'
    return f'fk = FK_INT; fw = {w}; fa = (uint64_t)({res}); fb = 0;'


def clab(l):
    return 'L' + l.replace('.', '_')


lines_out = []
labels_used = set()
for line in body:
    s = line.split('#')[0].rstrip()
    m = re.fullmatch(r'\s*([._A-Za-z0-9]+):\s*', s)
    if m:
        lines_out.append(('label', m.group(1)))
        continue
    if not s.strip() or s.strip().startswith('.'):
        continue
    t = s.strip()
    mn, _, rest = t.partition(' ')
    ops = split_ops(rest.strip())
    c = None
    if mn == 'endbr64':
        continue
    if mn in ('movsd',):
        a, b = ops
        if is_xmm(a) and is_xmm(b):
            c = f'{x(b)} = {x(a)};'
        elif is_xmm(b):
            c = f'{x(b)} = ldd({addr(a)});'
        else:
            c = f'std_({addr(b)}, {x(a)});'
    elif mn == 'movaps':
        c = f'{x(ops[1])} = {x(ops[0])};'
    elif mn in ('pxor', 'xorps') and ops[0] == ops[1]:
        c = f'{x(ops[1])} = fzero;'
    elif mn == 'xorps':
        a, b = ops
        c = f'{x(b)} = bitsd(dbits({x(b)}) ^ ld64({addr(a)}));'
    elif mn in ('addsd', 'subsd', 'mulsd', 'divsd'):
        a, b = ops
        o = {'add': '+', 'sub': '-', 'mul': '*', 'div': '/'}[mn[:3]]
        src_ = x(a) if is_xmm(a) else f'ldd({addr(a)})'
        c = f'{x(b)} = {x(b)} {o} {src_};'
    elif mn in ('comisd', 'ucomisd'):
        a, b = ops
        src_ = x(a) if is_xmm(a) else f'ldd({addr(a)})'
        c = f'fk = {"FK_COMI" if mn == "comisd" else "FK_UCOMI"}; fx = {x(b)}; fy = {src_};'
    elif mn in ('movq', 'movl', 'movw', 'movb'):
        w = SFX[mn[-1]]
        a, b = ops
        c = wr(b, w, rd(a, w))
    elif mn == 'movzwl':
        c = wr(ops[1], 32, rd(ops[0], 16))
    elif mn == 'movslq':
        c = wr(ops[1], 64, f'(int64_t)({ST[32]}){rd(ops[0], 32)}')
    elif mn == 'lea':
        r, _ = ALIAS[ops[1][1:]]
        c = f'{r} = {addr(ops[0])};'
    elif re.fullmatch(r'(and|or|xor|add|sub)[lbq]', mn):
        w = SFX[mn[-1]]
        o = {'and': '&', 'or': '|', 'xor': '^', 'add': '+', 'sub': '-'}[mn[:-1]]
        a, b = ops
        if mn[:-1] == 'xor' and a == b and not is_mem(a):
            res = '0'
        else:
            res = f'({UT[w]})({rd(b, w)} {o} {rd(a, w)})'
        if mn[:-1] == 'sub':
            fl = flags_int(None, w, rd(b, w), rd(a, w))
            c = f'{fl} ' + wr(b, w, res)
        else:
            c = f'{{ {UT[w]} t_ = {res}; {wr(b, w, "t_")} {flags_int("t_", w)} }}'
    elif re.fullmatch(r'(shl|shr|sar)[lbq]', mn):
        w = SFX[mn[-1]]
        a, b = ops
        n = int(a[1:], 0)
        if mn.startswith('sar'):
            res = f'({UT[w]})(({ST[w]}){rd(b, w)} >> {n})'
        else:
            res = f'({UT[w]})({rd(b, w)} {"<<" if mn.startswith("shl") else ">>"} {n})'
        c = f'{{ {UT[w]} t_ = {res}; {wr(b, w, "t_")} fk = FK_BAD; }}'
    elif mn == 'negl':
        c = f'{{ uint32_t t_ = (uint32_t)(0u - {rd(ops[0], 32)}); {wr(ops[0], 32, "t_")} fk = FK_BAD; }}'
    elif re.fullmatch(r'(test|cmp)[lbq]', mn):
        w = SFX[mn[-1]]
        a, b = ops
        if mn.startswith('test'):
            c = f'{{ {UT[w]} t_ = ({UT[w]})({rd(b, w)} & {rd(a, w)}); {flags_int("t_", w)} }}'
        else:
            c = flags_int(None, w, rd(b, w), rd(a, w))
    elif mn == 'sete':
        c = wr(ops[0], 8, 'cc_e(fk, fw, fa, fb, fx, fy)')
    elif mn == 'cmovl':
        c = wr(ops[1], 32, f'cc_l4(fk, fw, fa, fb) ? {rd(ops[0], 32)} : {rd(ops[1], 32)}')
    elif mn == 'jmp':
        labels_used.add(ops[0]); c = f'goto {clab(ops[0])};'
    elif mn in ('je', 'jne', 'jb', 'jbe', 'ja', 'jae', 'jl', 'jle', 'jg', 'jge', 'jp', 'jnp'):
        labels_used.add(ops[0])
        c = f'if (cc_{mn[1:]}(fk, fw, fa, fb, fx, fy)) goto {clab(ops[0])};'
    elif mn == 'ret':
        c = 'return (int)(uint32_t)rax;'
    else:
        raise SystemExit(f'unsupported instruction: {t}')
    lines_out.append(('insn', c, t))

# --- output ---
print(f'/* Generated by functions/asin/tools/asm2c.py from {os.path.basename(asm_path)} ({sym}); do not edit. */')
for s_ in tables:
    words = extract(asm_path, s_)
    print(f'/* {s_} ({len(words) * 8} bytes); writable so that operations on it are not constant-folded */')
    print(f'HWYMATH_HIDDEN uint64_t {tab_name(s_)}_{cname}[{len(words)}] = {{')
    for i in range(0, len(words), 4):
        print('    ' + ' '.join(f'0x{w:016x}ULL,' for w in words[i:i + 4]))
    print('};')
    print(f'#define {tab_name(s_)} {tab_name(s_)}_{cname}')
params = ', '.join(f'{"const double" if i < nargs - 1 else "double"} *a{i}' for i in range(nargs))
print(f'\nint\n{cname}({params})\n{{')
print('    uint64_t rax = 0, rbx = 0, rcx = 0, rdx = 0, rsi = 0, rdi = 0, rbp = 0, r8 = 0, r9 = 0, r10 = 0, '
      'r11 = 0, r12 = 0, r13 = 0, r14 = 0, r15 = 0;')
print('    double ' + ', '.join(f'x{i} = 0' for i in range(16)) + ';')
print('    int fk = FK_BAD, fw = 32;')
print('    uint64_t fa = 0, fb = 0;')
print('    double fx = 0, fy = 0;')
print('    _Alignas(16) unsigned char stack[256];')
print('    uint64_t rsp = (uint64_t)(uintptr_t)(stack + 192);  /* the routine uses the red zone below %rsp */')
regs = ['rdi', 'rsi', 'rdx']
for i in range(nargs):
    print(f'    {regs[i]} = (uint64_t)(uintptr_t)a{i};')
print('    (void)rbx; (void)rbp; (void)r12; (void)r13; (void)r14; (void)r15; (void)rsp;')
used = {clab(l) for l in labels_used}
for item in lines_out:
    if item[0] == 'label':
        lab = clab(item[1])
        if lab in used:
            print(f'{lab}:')
    else:
        print(f'    {item[1]:<96} /* {re.sub(r"\s+", " ", item[2])} */')
print('}')
for s_ in tables:
    print(f'#undef {tab_name(s_)}')
