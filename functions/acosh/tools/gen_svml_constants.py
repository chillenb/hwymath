"""Print C declarations for the acosh_d_ha constants, taken bit-exactly from the SVML assembly.
The data layout is that of asinh (functions/asinh/tools/gen_svml_constants.py does the work)."""
import os, runpy, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..', '..', '..')
asm = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'third_party', 'SVML', 'linux', 'avx512', 'svml_z0_acosh_d_ha.s')
sys.argv = [sys.argv[0], asm, 'acosh']
runpy.run_path(os.path.join(HERE, '..', '..', 'asinh', 'tools', 'gen_svml_constants.py'), run_name='__main__')
