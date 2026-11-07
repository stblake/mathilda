# Holdout-A D9 family (experiment 214) references, method A: EXACT rational
# propagation of the start distribution through `steps` steps with Python
# fractions over a SPARSE neighbour list (no matrix).  Method B: hd9_ref.wl
# (exact MatrixPower of the SparseArray transition matrix); hd9.7 also has the
# closed form 2^-8 Sum_k C(8,k) (-1)^k (1 - k/4)^30 (Krawtchouk / Fourier on Z_2^8).
from fractions import Fraction as Fr
from math import comb
import mpmath as mp
mp.mp.dps = 80

def tri(n):
    pts = [(r, c) for r in range(n) for c in range(r + 1)]
    nb = lambda p: [q for q in ((p[0], p[1]-1), (p[0], p[1]+1), (p[0]-1, p[1]-1), (p[0]-1, p[1]),
                                (p[0]+1, p[1]), (p[0]+1, p[1]+1)) if 0 <= q[1] <= q[0] < n]
    return pts, nb
def grid(n, king=False):
    pts = [(i, j) for i in range(1, n+1) for j in range(1, n+1)]
    ds = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0) and (king or a*b == 0)]
    nb = lambda p: [(p[0]+a, p[1]+b) for a, b in ds if 1 <= p[0]+a <= n and 1 <= p[1]+b <= n]
    return pts, nb
def cube(n):
    pts = [(i, j, k) for i in range(1, n+1) for j in range(1, n+1) for k in range(1, n+1)]
    ds = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    nb = lambda p: [tuple(a+b for a, b in zip(p, d)) for d in ds if all(1 <= a+b <= n for a, b in zip(p, d))]
    return pts, nb
def hyper(d):
    return list(range(2 ** d)), (lambda p: [p ^ (1 << i) for i in range(d)])

def walk(pts, nb, start, steps, target, lazy=False, absorb=lambda p: False):
    T = {p: ([(p, Fr(1))] if absorb(p) else
             ([(p, Fr(1, 2))] + [(q, Fr(1, 2*len(nb(p)))) for q in nb(p)] if lazy else
              [(q, Fr(1, len(nb(p)))) for q in nb(p)])) for p in pts}
    v = {start: Fr(1)}
    for _ in range(steps):
        w = {}
        for p, a in v.items():
            for q, b in T[p]:
                w[q] = w.get(q, 0) + a * b
        v = w
    return sum(a for p, a in v.items() if target(p))

R = {}
P, N = tri(10);  R[1] = walk(P, N, (0, 0), 20, lambda p: p[0] == 9)
P, N = tri(30);  R[2] = walk(P, N, (0, 0), 60, lambda p: p == (0, 0))
P, N = tri(20);  R[3] = walk(P, N, (12, 6), 50, lambda p: p[1] == 0 or p[1] == p[0] or p[0] == 19)
P, N = grid(20); R[4] = walk(P, N, (1, 1), 60, lambda p: p[0] == 20)
P, N = grid(12, king=True); R[5] = walk(P, N, (1, 1), 30, lambda p: p[0] == p[1])
P, N = cube(6);  R[6] = walk(P, N, (1, 1, 1), 41, lambda p: p == (6, 6, 6))
P, N = hyper(8); R[7] = walk(P, N, 0, 30, lambda p: p == 255)
P, N = tri(30);  R[8] = walk(P, N, (0, 0), 60, lambda p: p[0] == 29, lazy=True)
P, N = tri(30);  R[9] = walk(P, N, (0, 0), 60, lambda p: p[0] == 29, absorb=lambda p: p[0] == 29)
P, N = tri(40);  R[10] = walk(P, N, (0, 0), 80, lambda p: p[0] == 39)
c7 = Fr(sum(comb(8, k) * (-1) ** k * Fr(4 - k, 4) ** 30 for k in range(9)), 2 ** 8)
print("hd9.7 closed form equal:", c7 == R[7])
for k in range(1, 11):
    q = R[k]
    print("hd9.%d  %s  num %d digits / den %d digits" % (k, mp.nstr(mp.mpf(q.numerator) / q.denominator, 70),
                                                         len(str(q.numerator)), len(str(q.denominator))))
    print("   exact %d/%d" % (q.numerator, q.denominator))
