"""Reference for d9 (triangular walk): exact rational by Python fractions,
pushing the probability vector 60 steps on the 465-point lattice."""
from fractions import Fraction
from mpmath import mp, mpf
n = 30
idx = {}
for r in range(n):
    for c in range(r + 1):
        idx[(r, c)] = len(idx)
def nbrs(r, c):
    for (a, b) in ((r, c-1), (r, c+1), (r-1, c-1), (r-1, c), (r+1, c), (r+1, c+1)):
        if 0 <= b <= a < n:
            yield (a, b)
v = {(0, 0): Fraction(1)}
for _ in range(60):
    w = {}
    for p, pr in v.items():
        nb = list(nbrs(*p))
        q = pr / len(nb)
        for t in nb:
            w[t] = w.get(t, 0) + q
    v = w
p = sum(pr for (r, c), pr in v.items() if r == n - 1)
print("exact", p.numerator, "/", p.denominator)
mp.dps = 70
print(mp.nstr(mpf(p.numerator) / p.denominator, 65))
