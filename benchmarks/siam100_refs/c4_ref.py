"""c4 reference (method 1): Newton on grad f = 0 (mpmath.findroot, analytic
gradient), dps 90, started from the global minimiser located by
scipy differential_evolution; f evaluated at the critical point."""
import mpmath as mp
from scipy.optimize import differential_evolution
import math
mp.mp.dps = 90
def fm(x, y): return (math.exp(math.sin(50*x)) + math.sin(60*math.exp(y)) + math.sin(70*math.sin(x))
                      + math.sin(math.sin(80*y)) - math.sin(10*(x+y)) + (x*x+y*y)/4)
r = differential_evolution(lambda v: fm(*v), [(-1, 1), (-1, 1)], tol=1e-14, seed=1, maxiter=3000, popsize=60)
f = lambda x, y: (mp.exp(mp.sin(50*x)) + mp.sin(60*mp.exp(y)) + mp.sin(70*mp.sin(x))
                  + mp.sin(mp.sin(80*y)) - mp.sin(10*(x+y)) + (x*x+y*y)/4)
gx = lambda x, y: 50*mp.cos(50*x)*mp.exp(mp.sin(50*x)) + 70*mp.cos(x)*mp.cos(70*mp.sin(x)) - 10*mp.cos(10*(x+y)) + x/2
gy = lambda x, y: 60*mp.exp(y)*mp.cos(60*mp.exp(y)) + 80*mp.cos(80*y)*mp.cos(mp.sin(80*y)) - 10*mp.cos(10*(x+y)) + y/2
x, y = mp.findroot([gx, gy], (mp.mpf(r.x[0]), mp.mpf(r.x[1])))
print(r.x, r.fun)
print(mp.nstr(x, 70), mp.nstr(y, 70))
print(mp.nstr(f(x, y), 70))
