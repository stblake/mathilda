# Writes d3_start.wl: the double-precision global minimiser (scipy BFGS with the
# analytic gradient from the Fibonacci start, as in d3_ref.py) for d3_ref.wl.
import numpy as np, scipy.optimize as so
src = open('d3_ref.py').read().split('mp.mp.dps = 130')[0].split('rng = np.random')[0]
exec(src)
with open('d3_start.wl', 'w') as fh:
    fh.write('d3start = {' + ', '.join(repr(float(t)) for t in r.x) + '};\n')
print(r.fun)
