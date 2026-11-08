"""Reference A for holdout-B family D18 (exp 253): smallest Dirichlet
eigenvalue of -Laplacian on domains with a closed-form lambda_1.

  bd18.1  unit square                     2 Pi^2
  bd18.2  right isosceles triangle, legs 1   5 Pi^2   (mode sin(2pi x)sin(pi y) - sin(pi x)sin(2pi y))
  bd18.3  hemi-equilateral (30-60-90) triangle = half of the equilateral
          triangle of side 1: the equilateral's first antisymmetric mode,
          (16 Pi^2/9)(m^2 + m n + n^2) at (m,n) = (1,2): 112 Pi^2/9
  bd18.4  unit disk                        j_{0,1}^2
  bd18.5  sector of angle 3Pi/2, radius 1  j_{2/3,1}^2   (order Pi/angle)
  bd18.6  annulus 1/2 < r < 1              k^2, k the first root of
                                           J0(k/2) Y0(k) - J0(k) Y0(k/2)
  bd18.7  unit cube                        3 Pi^2
  bd18.8  unit ball                        Pi^2          (j_{1/2,1} = Pi)
  bd18.9  cylinder radius 1, height 1      j_{0,1}^2 + Pi^2
  bd18.10 sector of angle 7Pi/4, radius 1  j_{4/7,1}^2
Method A (this file): mpmath besseljzero / findroot at 80 digits.
Method B: bd18_ref.gp (PARI solve() on besselj / besseln, independent code).
Plus a low-precision check of each closed-form CLAIM for the polygons
(bd18.1-3): P1 finite elements on uniformly refined meshes + Richardson
extrapolation (scipy), which must agree to ~6 digits.
Usage: python3.11 bd18_ref.py
"""
import numpy as np
import mpmath as mp
import scipy.sparse as sp
import scipy.sparse.linalg as sla

mp.mp.dps = 80
pi = mp.pi
j01 = mp.besseljzero(0, 1)
ann = mp.findroot(lambda k: mp.besselj(0, k / 2) * mp.bessely(0, k) - mp.besselj(0, k) * mp.bessely(0, k / 2), 6.2)
vals = {
    "bd18.1": 2 * pi ** 2, "bd18.2": 5 * pi ** 2, "bd18.3": 112 * pi ** 2 / 9,
    "bd18.4": j01 ** 2, "bd18.5": mp.besseljzero(mp.mpf(2) / 3, 1) ** 2, "bd18.6": ann ** 2,
    "bd18.7": 3 * pi ** 2, "bd18.8": pi ** 2, "bd18.9": j01 ** 2 + pi ** 2,
    "bd18.10": mp.besseljzero(mp.mpf(4) / 7, 1) ** 2}
for k in sorted(vals, key=lambda s: int(s.split(".")[1])):
    print(k, mp.nstr(vals[k], 70), flush=True)


# ---- low-precision FEM check of the polygon closed forms ----
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


def lam1(P, T):
    P = np.array(P); n = len(P)
    I, J, K, M = [], [], [], []
    for t in T:
        x = P[list(t)]
        B = np.array([x[1] - x[0], x[2] - x[0]]).T
        area = abs(np.linalg.det(B)) / 2
        G = np.linalg.inv(B).T @ np.array([[-1, 1, 0], [-1, 0, 1]])
        Ke = area * G.T @ G
        Me = area / 12 * (np.ones((3, 3)) + np.eye(3))
        for r in range(3):
            for s in range(3):
                I.append(t[r]); J.append(t[s]); K.append(Ke[r, s]); M.append(Me[r, s])
    Kg = sp.csr_matrix((K, (I, J)), shape=(n, n)); Mg = sp.csr_matrix((M, (I, J)), shape=(n, n))
    # boundary nodes: edges used by exactly one triangle
    from collections import Counter
    ec = Counter(tuple(sorted(e)) for t in T for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])))
    bnd = {v for e, c in ec.items() if c == 1 for v in e}
    free = np.array([i for i in range(n) if i not in bnd])
    Kf = Kg[free][:, free]; Mf = Mg[free][:, free]
    return sla.eigsh(Kf, k=1, M=Mf, sigma=0, which="LM")[0][0]


s3 = np.sqrt(3) / 2
for lab, P, T in [("bd18.1", [(0, 0), (1, 0), (1, 1), (0, 1)], [(0, 1, 2), (0, 2, 3)]),
                  ("bd18.2", [(0, 0), (1, 0), (0, 1)], [(0, 1, 2)]),
                  ("bd18.3", [(0, 0), (0.5, 0), (0.5, s3)], [(0, 1, 2)])]:
    ls = []
    for lev in range(7):
        if lev >= 4:
            ls.append(lam1(P, T))
        if lev < 6:
            P, T = refine(P, T)
    r1 = [(4 * ls[i + 1] - ls[i]) / 3 for i in range(len(ls) - 1)]
    r2 = (16 * r1[1] - r1[0]) / 15
    print("FEM check", lab, "levels 4..6:", ["%.8f" % v for v in ls], "Richardson %.9f" % r2,
          "closed form %.9f" % float(vals[lab]), flush=True)
