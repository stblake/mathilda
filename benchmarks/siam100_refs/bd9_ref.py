# Holdout B, family D9 -- method A: exact rational propagation of the
# start distribution (sparse dict of Fractions), then the target mass.
# Method B: bd9_ref.gp builds the graphs itself (independent code) and
# propagates exact rationals in PARI/gp; closed forms for bd9.7 and bd9.8.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd9_common import walks, trans
from fractions import Fraction as Fr
from math import comb
import mpmath as mp
mp.mp.dps = 70
W = walks()
for lab, (nb, T, s, tg, n) in W.items():
    P = trans(nb, T)
    v = {s: Fr(1)}
    for _ in range(n):
        w = {}
        for i, p in v.items():
            for j, q in P[i]:
                w[j] = w.get(j, 0) + p * q
        v = w
    x = sum(v.get(j, 0) for j in tg)
    print(lab, mp.nstr(mp.mpf(x.numerator) / x.denominator, 66), "exact %d/%d" % (x.numerator, x.denominator) if len(str(x.denominator)) < 400 else "")

# closed forms
f7 = sum((-1) ** k * comb(10, k) * Fr(10 - 2 * k, 10) ** 40 for k in range(11)) / 1024
q7 = mp.mpf(f7.numerator) / f7.denominator
print("bd9.7 closed", mp.nstr(q7, 66))
q8 = mp.mpf(comb(60, 30) + 2 * comb(60, 47)) / 2 ** 60
print("bd9.8 closed", mp.nstr(q8, 66))
