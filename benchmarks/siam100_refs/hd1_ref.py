"""Holdout A, family hd1 (experiment 206): Kempner-type sums over digit-restricted n.

Generic Schmelzer-Baillie power-sum recursion over a finite digit automaton:
  psi_t(k, q) = sum of n^-t over admissible k-digit n ending in automaton state q,
  (b n + d)^-t = sum_w binom(-t, w) d^w b^(-t-w) n^(-t-w),
so psi(k+1) = M psi(k) with a constant matrix M (exponents t = s .. s+J).
Explicit DFS enumeration of the admissible n with <= K0 digits, then the tail
  sum_{k>K0} psi_s(k, .) = e^T (I - M)^-1 M psi(K0)     (closed form).
Usage: hd1_ref.py K0 J dps [problem ...]"""
import sys
from mpmath import mp, mpf, matrix, binomial, lu_solve, nstr


def kmp(pattern, b):
    """Substring-avoidance automaton: state = length of the matched prefix."""
    P = [int(c) for c in pattern]
    m = len(P)

    def delta(q, d):
        s = P[:q] + [d]
        for L in range(min(len(s), m), -1, -1):
            if L == 0 or s[len(s) - L:] == P[:L]:
                return None if L == m else L
    return list(range(m)), 0, delta


def nodigit(dset):
    return [0], 0, (lambda q, d: None if d in dset else 0)


def nondecr():
    # state = last digit (0 = start)
    return list(range(10)), 0, (lambda q, d: d if d >= max(q, 1) else None)


PROBLEMS = {
    1: (10, nodigit({9}), 1),
    2: (10, nodigit({0}), 1),
    3: (10, kmp("00", 10), 1),
    4: (10, kmp("314", 10), 1),
    5: (2, kmp("11", 2), 1),
    6: (3, nodigit({2}), 1),
    7: (10, nodigit({9}), 2),
    8: (10, nodigit({0, 2, 4, 6, 8}), 1),
    9: (10, nondecr(), 1),
    10: (10, kmp("2026", 10), 1),
}


def kempner(b, aut, s, K0, J):
    states, q0, delta = aut
    ns = len(states)
    sidx = {q: i for i, q in enumerate(states)}
    T = [mpf(s + j) for j in range(J + 1)]
    psi = [[mpf(0)] * (J + 1) for _ in range(ns)]
    total = mpf(0)
    # DFS over admissible numbers with <= K0 digits
    stack = []
    for d in range(1, b):
        q = delta(q0, d)
        if q is not None:
            stack.append((d, q, 1))
    while stack:
        n, q, k = stack.pop()
        x = mpf(1) / n
        total += x ** s
        if k == K0:
            p = x ** s
            for j in range(J + 1):
                psi[sidx[q]][j] += p
                p *= x
        else:
            for d in range(b):
                q2 = delta(q, d)
                if q2 is not None:
                    stack.append((n * b + d, q2, k + 1))
    N = ns * (J + 1)
    M = matrix(N, N)
    for q in states:
        for d in range(b):
            q2 = delta(q, d)
            if q2 is None:
                continue
            for j in range(J + 1):
                t = T[j]
                for w in range(0, J - j + 1):
                    c = binomial(-t, w) * mpf(d) ** w / mpf(b) ** (t + w)
                    if c != 0:
                        M[sidx[q2] * (J + 1) + j, sidx[q] * (J + 1) + j + w] += c
    v = matrix([psi[i][j] for i in range(ns) for j in range(J + 1)])
    I = matrix(N, N)
    for i in range(N):
        I[i, i] = 1
    y = lu_solve(I - M, M * v)
    for i in range(ns):
        total += y[i * (J + 1)]
    return total


if __name__ == "__main__":
    K0, J, dps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    mp.dps = dps
    which = [int(a) for a in sys.argv[4:]] or sorted(PROBLEMS)
    for k in which:
        b, aut, s = PROBLEMS[k]
        K = K0 if b == 10 else {2: 4 * K0, 3: 2 * K0 + 1}[b]
        print(k, nstr(kempner(b, aut, s, K, J), dps - 10), flush=True)
