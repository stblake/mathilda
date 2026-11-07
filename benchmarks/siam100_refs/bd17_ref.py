"""Reference A for holdout-B family D17 (exp 252): length of the nodal set
(zero set in the open domain) of the k-th Dirichlet eigenfunction, on domains
where that eigenvalue is SIMPLE and the eigenfunction separable.

  bd17.1  rectangle [0,2]x[0,1], k=2 ((m,n) = (2,1): sin(pi x) sin(pi y)): x=1: 1
  bd17.2  unit square, k=4 (8 pi^2, sin 2pi x sin 2pi y): x=1/2, y=1/2: 2
  bd17.3  unit disk, k=6 (j02^2, J0(j02 r)): circle r = j01/j02: 2 pi j01/j02
  bd17.4  upper half disk, k=4 (j12^2, J1(j12 r) sin t): semicircle r = j11/j12: pi j11/j12
  bd17.5  triangle (0,0),(1,0),(1,1), k=2 ((3,1): sin3a sin b - sin a sin3b
          = 4 sin a sin b (sin^2 b - sin^2 a)): segment x+y=1: sqrt2/2
  bd17.6  same triangle, k=3 ((3,2): = 2 sin a sin b (cos a - cos b)(4 cos a cos b + 1)):
          the curve cos(pi x) cos(pi y) = -1/4 from (arccos(-1/4)/pi, 0) to (1, arccos(1/4)/pi)
  bd17.7  quarter disk, k=3 (j22^2, J2(j22 r) sin 2t): quarter circle r = j21/j22: (pi/2) j21/j22
  bd17.8  [0,3/2]x[0,1], u=0 except NEUMANN on y=1, k=3 (sin(2pi x/3) sin(3pi y/2)): y=2/3: 3/2
  bd17.9  equilateral triangle of side 1, k=4 ((2,2), 64 pi^2/3): the medial triangle: 3/2
  bd17.10 unit disk, k=15 (j03^2, J0(j03 r)): circles j01/j03, j02/j03: 2 pi (j01+j02)/j03
Method A (this file): mpmath at 80 digits; bd17.6 by tanh-sinh quadrature of
the two halves of the curve, one parametrised by y, the other by x.
Method B: bd17_ref.gp (PARI: Bessel zeros by solve(); bd17.6 as TWICE the
half-curve x in [2/3, 1] -- the curve is symmetric under (x,y) -> (1-y,1-x) --
by intnum).  Index / simplicity claims for the polygons are checked by P1 FEM
(first 6 eigenvalues, scipy) against pi^2 (m^2+n^2) resp. (16 pi^2/9)(m^2+mn+n^2).
Usage: python3.11 bd17_ref.py
"""
import numpy as np
import mpmath as mp
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from collections import Counter

mp.mp.dps = 80
pi = mp.pi
z = lambda m, k: mp.besseljzero(m, k)


def yx(x):
    return mp.acos(-1 / (4 * mp.cos(pi * x))) / pi


def dy(x):
    u = -1 / (4 * mp.cos(pi * x))
    return mp.sin(pi * x) / (4 * mp.cos(pi * x) ** 2 * mp.sqrt(1 - u ** 2))


x0 = mp.acos(-mp.mpf(1) / 4) / pi
# the curve meets the anti-diagonal at (2/3, 1/3).  Piece 1, from (x0, 0) to (2/3, 1/3), is
# parametrised by y (x(y) = arccos(-1/(4 cos pi y))/pi is smooth there; y(x) has a square-root
# endpoint singularity at x0, which costs a tanh-sinh rule ~35 digits through cancellation in
# 1 - u^2); piece 2, from (2/3, 1/3) to (1, y1), by x.  By the symmetry (x,y) -> (1-y,1-x) the two
# are equal -- printed as a check; method B (gp) uses 2 x piece 2 only.
L6a = mp.quad(lambda y: mp.sqrt(1 + dy(y) ** 2), [0, mp.mpf(1) / 3])   # x'(y) has the same formula
L6b = mp.quad(lambda x: mp.sqrt(1 + dy(x) ** 2), [mp.mpf(2) / 3, 1])
print("bd17.6 pieces", mp.nstr(L6a, 30), mp.nstr(L6b, 30))
L6 = L6a + L6b
assert abs(mp.cos(pi * 1) * mp.cos(pi * yx(mp.mpf(1)))+ mp.mpf(1)/4) < mp.mpf(10)**-70
vals = {"bd17.1": mp.mpf(1), "bd17.2": mp.mpf(2), "bd17.3": 2 * pi * z(0, 1) / z(0, 2),
        "bd17.4": pi * z(1, 1) / z(1, 2), "bd17.5": mp.sqrt(2) / 2, "bd17.6": L6,
        "bd17.7": pi / 2 * z(2, 1) / z(2, 2), "bd17.8": mp.mpf(3) / 2, "bd17.9": mp.mpf(3) / 2,
        "bd17.10": 2 * pi * (z(0, 1) + z(0, 2)) / z(0, 3)}
for k in sorted(vals, key=lambda s: int(s.split(".")[1])):
    print(k, mp.nstr(vals[k], 70), flush=True)

# disk ordering (bd17.3, bd17.10): zeros j_{m,k} with multiplicity 2 for m >= 1
zs = sorted([(float(z(m, k)), m, k) for m in range(0, 8) for k in range(1, 5)])
cnt = 0
for v, m, k in zs:
    cnt += 1 if m == 0 else 2
    if v < 9:
        print("disk: j_{%d,%d} = %.4f  eigenvalue index up to %d" % (m, k, v, cnt))
# half disk J_m sin(m t), m >= 1 (simple); quarter disk J_{2m} sin(2 m t)
print("half disk order:", sorted([(round(float(z(m, k)), 4), m, k) for m in range(1, 6) for k in range(1, 3)])[:5])
print("quarter disk order:", sorted([(round(float(z(2 * m, k)), 4), 2 * m, k) for m in range(1, 4) for k in range(1, 3)])[:4])


# ---- P1 FEM check of the polygon spectra ----
def refine(P, T):
    P = list(map(tuple, P)); idx = {p: i for i, p in enumerate(P)}; T2 = []
    def mid(a, b):
        m = ((P[a][0] + P[b][0]) / 2, (P[a][1] + P[b][1]) / 2)
        if m not in idx:
            idx[m] = len(P); P.append(m)
        return idx[m]
    for a, b, c in T:
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        T2 += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
    return P, T2


def eigs(P, T, k, neumann=lambda p: False):
    P = np.array(P); n = len(P); I, J, K, M = [], [], [], []
    for t in T:
        x = P[list(t)]
        B = np.array([x[1] - x[0], x[2] - x[0]]).T
        area = abs(np.linalg.det(B)) / 2
        G = np.linalg.inv(B).T @ np.array([[-1, 1, 0], [-1, 0, 1]])
        Ke = area * G.T @ G; Me = area / 12 * (np.ones((3, 3)) + np.eye(3))
        for r in range(3):
            for s in range(3):
                I.append(t[r]); J.append(t[s]); K.append(Ke[r, s]); M.append(Me[r, s])
    Kg = sp.csr_matrix((K, (I, J)), shape=(n, n)); Mg = sp.csr_matrix((M, (I, J)), shape=(n, n))
    ec = Counter(tuple(sorted(e)) for t in T for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])))
    bnd = set()
    for e, c in ec.items():
        if c == 1 and not (neumann(P[e[0]]) and neumann(P[e[1]])):
            bnd |= set(e)
    free = np.array([i for i in range(n) if i not in bnd])
    return np.sort(sla.eigsh(Kg[free][:, free], k=k, M=Mg[free][:, free], sigma=0, which="LM")[0])


s3 = np.sqrt(3) / 2
P2 = float(pi) ** 2
for lab, P, T, nm, ref in [
        ("rect 2x1 (bd17.1)", [(0, 0), (2, 0), (2, 1), (0, 1)], [(0, 1, 2), (0, 2, 3)], lambda p: False,
         sorted(P2 * (m * m / 4 + n * n) for m in range(1, 8) for n in range(1, 4))[:6]),
        ("square (bd17.2)", [(0, 0), (1, 0), (1, 1), (0, 1)], [(0, 1, 2), (0, 2, 3)], lambda p: False,
         sorted(P2 * (m * m + n * n) for m in range(1, 5) for n in range(1, 5))[:6]),
        ("triangle (bd17.5/6)", [(0, 0), (1, 0), (1, 1)], [(0, 1, 2)], lambda p: False,
         sorted(P2 * (m * m + n * n) for m in range(2, 7) for n in range(1, m))[:6]),
        ("mixed rect (bd17.8)", [(0, 0), (1.5, 0), (1.5, 1), (0, 1)], [(0, 1, 2), (0, 2, 3)],
         lambda p: abs(p[1] - 1) < 1e-12, sorted(P2 * ((2 * m / 3) ** 2 + (n - 0.5) ** 2) for m in range(1, 6) for n in range(1, 4))[:6]),
        ("equilateral (bd17.9)", [(0, 0), (1, 0), (0.5, s3)], [(0, 1, 2)], lambda p: False,
         sorted(16 * P2 / 9 * (m * m + m * n + n * n) for m in range(1, 6) for n in range(1, 6))[:6])]:
    for _ in range(6):
        P, T = refine(P, T)
    print("FEM", lab, ["%.4f" % v for v in eigs(P, T, 6, nm)], "\n   predicted", ["%.4f" % v for v in ref], flush=True)
