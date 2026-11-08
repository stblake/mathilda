# Holdout B, family D9: the ten walks (graph, start, steps, target set).
# Each walk moves from its current vertex to a neighbour chosen uniformly
# (bd9.9: the gambler's-ruin chain with absorbing ends).  Vertices are
# numbered 0.. in the order listed here, which is the order every .m/.py/.mac
# /.input half uses (1-based there).
from fractions import Fraction as Fr

def tri(n):
    pts = [(r, c) for r in range(n) for c in range(r + 1)]
    ix = {p: i for i, p in enumerate(pts)}
    nb = []
    for (r, c) in pts:
        l = [(a, b) for (a, b) in ((r, c-1), (r, c+1), (r-1, c-1), (r-1, c), (r+1, c), (r+1, c+1)) if 0 <= b <= a < n]
        nb.append([ix[p] for p in l])
    return pts, nb, ix

def grid(dims, moves):
    import itertools
    pts = list(itertools.product(*[range(d) for d in dims]))
    ix = {p: i for i, p in enumerate(pts)}
    nb = []
    for p in pts:
        l = [tuple(a + b for a, b in zip(p, m)) for m in moves]
        nb.append([ix[q] for q in l if q in ix])
    return pts, nb, ix

def walks():
    W = {}
    # 1: path 0..20, start 0, 40 steps, end in 11..20
    pts, nb, ix = grid([21], [(1,), (-1,)])
    W["bd9.1"] = (nb, None, ix[(0,)], [ix[(i,)] for i in range(11, 21)], 40)
    # 2: 10x10 grid (4 neighbours), start corner (0,0), 30 steps, end in row 9
    pts, nb, ix = grid([10, 10], [(1, 0), (-1, 0), (0, 1), (0, -1)])
    W["bd9.2"] = (nb, None, ix[(0, 0)], [ix[(9, j)] for j in range(10)], 30)
    # 3: triangular lattice, 20 per side, start top, 40 steps, end in bottom row
    pts, nb, ix = tri(20)
    W["bd9.3"] = (nb, None, ix[(0, 0)], [ix[(19, c)] for c in range(20)], 40)
    # 4: 6x6x6 cube grid (6 neighbours), start corner, 24 steps, end on face x = 5
    mv = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
    pts, nb, ix = grid([6, 6, 6], mv)
    W["bd9.4"] = (nb, None, ix[(0, 0, 0)], [ix[p] for p in pts if p[0] == 5], 24)
    # 5: knight on 8x8 from a1 = (0,0), 20 moves, end on h8 = (7,7)
    mv = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
    pts, nb, ix = grid([8, 8], mv)
    W["bd9.5"] = (nb, None, ix[(0, 0)], [ix[(7, 7)]], 20)
    # 6: king on 8x8 from a1, 30 moves, end on rank 8 (second coordinate 7)
    mv = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
    pts, nb, ix = grid([8, 8], mv)
    W["bd9.6"] = (nb, None, ix[(0, 0)], [ix[(i, 7)] for i in range(8)], 30)
    # 7: hypercube Q_10 (vertex v ~ v xor 2^i), start 0, 40 steps, end at 1023
    nb = [[v ^ (1 << i) for i in range(10)] for v in range(1024)]
    W["bd9.7"] = (nb, None, 0, [1023], 40)
    # 8: cycle C_17, start 0, 60 steps, end at 0
    nb = [[(v + 1) % 17, (v - 1) % 17] for v in range(17)]
    W["bd9.8"] = (nb, None, 0, [0], 60)
    # 9: gambler's ruin on 0..30, up 2/5, down 3/5, 0 and 30 absorbing;
    #    start 15, P(at 30 after 200 steps)
    T = []
    for v in range(31):
        if v in (0, 30):
            T.append([(v, Fr(1))])
        else:
            T.append([(v + 1, Fr(2, 5)), (v - 1, Fr(3, 5))])
    W["bd9.9"] = (None, T, 15, [30], 200)
    # 10: triangular lattice, 40 per side (820 points), 100 steps, bottom row
    pts, nb, ix = tri(40)
    W["bd9.10"] = (nb, None, ix[(0, 0)], [ix[(39, c)] for c in range(40)], 100)
    return W

def trans(nb, T):
    if T is not None:
        return T
    return [[(j, Fr(1, len(l))) for j in l] for l in nb]
