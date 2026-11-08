"""Holdout B, family C7 (experiment 232): the ten sparse matrices, built
exactly as the benchmark halves build them.  Each builder returns
(n, entries, b, idx): entries = dict {(i, j): int} (0-based), right-hand side
b (list of ints), and the solution component idx (0-based) that is asked for.
Shared by bc7_ref.py (references) and the checks."""
import numpy as np


def primes(n):
    N = 20 * n + 100
    s = np.ones(N, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0][:n]]


def banded(n, offs, sub=1, sup=1):
    P = primes(n)
    E = {(i, i): P[i] for i in range(n)}
    for o in offs:
        for i in range(n - o):
            E[(i, i + o)] = sup
            E[(i + o, i)] = sub
    return E


def unit(n, k=0):
    b = [0] * n; b[k] = 1; return b


def lap1(N):
    E = {}
    for i in range(N):
        E[(i, i)] = 2
        if i + 1 < N:
            E[(i, i + 1)] = -1; E[(i + 1, i)] = -1
    return E


def kronsum(Ts, N, shift):
    """sum_k I x .. x T x .. x I  + shift*I, T = lap1(N), d = len(Ts)."""
    d = Ts
    n = N ** d
    E = {}
    for idx in range(n):
        c = []; r = idx
        for _ in range(d):
            c.append(r % N); r //= N
        E[(idx, idx)] = 2 * d + shift
        st = 1
        for k in range(d):
            if c[k] + 1 < N:
                E[(idx, idx + st)] = -1; E[(idx + st, idx)] = -1
            st *= N
    return n, E


def build(k):
    if k == 1:   # tridiagonal, primes, n = 1000
        n = 1000; return n, banded(n, [1]), unit(n), 0
    if k == 2:   # 2-D 5-point Laplacian (100 x 100) + I, centre node (50,50)
        N = 100; n, E = kronsum(2, N, 1); return n, E, unit(n, 49 * N + 49), 49 * N + 49
    if k == 3:   # primes, 1 at |i-j| a Fibonacci number, n = 20000
        n = 20000; F = [1, 2]
        while F[-1] + F[-2] < n: F.append(F[-1] + F[-2])
        return n, banded(n, F), unit(n), 0
    if k == 4:   # primes, 1 at |i-j| a perfect square, n = 2000
        n = 2000; return n, banded(n, [j * j for j in range(1, 45) if j * j < n]), unit(n), 0
    if k == 5:   # primes, +1 above / -1 below at |i-j| = 2^k, n = 2000 (nonsymmetric)
        n = 2000; return n, banded(n, [2 ** j for j in range(11)], sub=-1, sup=1), unit(n), 0
    if k == 6:   # original structure, n = 5000, b = all ones, x_1
        n = 5000; return n, banded(n, [2 ** j for j in range(13) if 2 ** j < n]), [1] * n, 0
    if k == 7:   # Cayley graph of Z_p^* (p = 2003), generators 2, 3, 5; 7 on the diagonal
        p = 2003; E = {(i, i): 7 for i in range(p)}
        for v in range(1, p):
            for g in (2, 3, 5):
                w = g * v % p
                E[(v, w)] = 1; E[(w, v)] = 1
        return p, E, unit(p, 1), 1
    if k == 8:   # 3-D 7-point Laplacian on 30^3, corner node
        N = 30; n, E = kronsum(3, N, 0); return n, E, unit(n), 0
    if k == 9:   # 2-D Laplacian 150 x 150 minus 1/2 I (indefinite), centre (75,75)
        # returned SCALED BY 2 so that it is an integer matrix: 2A x = 2 e_c
        N = 150; n, E = kronsum(2, N, 0)
        E = {key: 2 * v for key, v in E.items()}
        for i in range(n):
            E[(i, i)] -= 1
        c = 74 * N + 74
        b = [0] * n; b[c] = 2
        return n, E, b, c
    if k == 10:  # complete binary tree, 2047 nodes: graph Laplacian + I, root
        n = 2047; E = {}
        for i in range(1, n + 1):
            deg = (1 if i > 1 else 0) + (2 if 2 * i <= n else 0)
            E[(i - 1, i - 1)] = deg + 1
            for c in (2 * i, 2 * i + 1):
                if c <= n:
                    E[(i - 1, c - 1)] = -1; E[(c - 1, i - 1)] = -1
        return n, E, unit(n), 0
