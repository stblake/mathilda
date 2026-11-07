"""Holdout B, family D1 reference, method A (mpmath): digit-restricted series
Sum_{n accepted} n^-s, for n accepted by a finite automaton reading the base-b
digits of n most-significant first (the Schmelzer-Baillie power-sum
recursion, generalised to any automaton and any s).

psi_q(k, j) = sum of n^-j over the admissible k-digit n that leave the
automaton in state q.  Appending a digit d to n:
   (b n + d)^-j = sum_w (-1)^w C(j+w-1, w) d^w b^(-j-w) n^(-j-w),
so psi(k+1, .) = M psi(k, .) with a CONSTANT matrix M (moments truncated at J).
Explicit enumeration for k <= K0, then the whole tail in closed form
   S = sum_{k<=K0} sum_{q accepting} psi_q(k, s) + e_acc^T (I - M)^-1 M psi(K0).

bd1.9 (Sum 2^-n over fibbinary n) is not a power sum: it is summed directly
to n = 2^9 here (error < 2^-500).

Usage: bd1_ref.py K0 J dps   -> prints one line per problem.
Independent method B: bd1_ref.gp (PARI, different K0 / J, power-iteration
cross-check where cheap, functional equation for bd1.9)."""
import sys
from mpmath import mp, mpf, matrix, binomial, lu_solve, nstr


def kmp_avoid(pat, base):
    """Automaton (states 0..len(pat)-1) of strings avoiding substring pat."""
    m = len(pat)

    def step(q, d):
        s = pat[:q] + d
        while s and not pat.startswith(s):
            s = s[1:]
        return len(s) if len(s) < m else None
    return step, list(range(m)), 0, lambda q: True


PROBLEMS = {}
# name: (base, digit chars, step(q,d)->q or None, states, q0, accept(q), s)
digs10 = "0123456789"
PROBLEMS["bd1.1"] = (10, 1, lambda q, d: None if d == "9" else 0, [0], 0, lambda q: True)
PROBLEMS["bd1.2"] = (10, 1, lambda q, d: None if d == "0" else 0, [0], 0, lambda q: True)
PROBLEMS["bd1.3"] = (10, 1, lambda q, d: 0 if d in "13579" else None, [0], 0, lambda q: True)
st, sts, q0, acc = kmp_avoid("314", 10)
PROBLEMS["bd1.4"] = (10, 1, st, sts, q0, acc)
PROBLEMS["bd1.5"] = (10, 1, lambda q, d: (q + (d == "9")) if q + (d == "9") <= 1 else None,
                     [0, 1], 0, lambda q: q == 1)
PROBLEMS["bd1.6"] = (10, 2, lambda q, d: None if d == "7" else 0, [0], 0, lambda q: True)
PROBLEMS["bd1.7"] = (2, 1, lambda q, d: None if (q == 1 and d == "1") else int(d), [0, 1], 0,
                     lambda q: True)
PROBLEMS["bd1.8"] = (3, 1, lambda q, d: None if d == "2" else 0, [0], 0, lambda q: True)
# bd1.10: no two equal adjacent decimal digits; state = last digit (10 = start)
PROBLEMS["bd1.10"] = (10, 1, lambda q, d: None if q == int(d) else int(d), list(range(10)), 10,
                      lambda q: True)


def todigits(n, b):
    s = ""
    while n:
        s = "0123456789"[n % b] + s
        n //= b
    return s


def run(base, s, step, states, q0, accept, K0, J):
    idx = {q: i for i, q in enumerate(states)}
    ns = len(states)
    psi = [[mpf(0)] * (J + 1) for _ in range(ns)]
    S = mpf(0)
    for n in range(1, base ** K0):
        ds = todigits(n, base)
        q = q0
        for ch in ds:
            q = step(q, ch)
            if q is None:
                break
        if q is None:
            continue
        x = mpf(1) / n
        if accept(q):
            S += x ** s
        if len(ds) == K0:
            p = x ** s
            for j in range(s, J + 1):
                psi[idx[q]][j] += p
                p *= x
    # moments j = s..J, index (state, j)
    nj = J - s + 1
    N = ns * nj
    pos = lambda q, j: idx[q] * nj + (j - s)
    M = matrix(N, N)
    for q in states:
        for dch in "0123456789"[:base]:
            l = step(q, dch)
            if l is None:
                continue
            d = int(dch)
            for j in range(s, J + 1):
                for w in range(0, J - j + 1):
                    c = (-1) ** w * binomial(j + w - 1, w) * mpf(d) ** w / mpf(base) ** (j + w)
                    M[pos(l, j), pos(q, j + w)] += c
    v = matrix([psi[idx[q]][j] for q in states for j in range(s, J + 1)])
    I = matrix(N, N)
    for i in range(N):
        I[i, i] = 1
    y = lu_solve(I - M, M * v)
    for q in states:
        if accept(q):
            S += y[pos(q, s)]
    return S


if __name__ == "__main__":
    K0, J, dps = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    mp.dps = dps
    only = sys.argv[4:] or None
    for name, (base, s, step, states, q0, accept) in PROBLEMS.items():
        if only and name not in only:
            continue
        k0 = K0 if base == 10 else {2: 4 * K0, 3: 2 * K0 + 1}[base]
        print(name, nstr(run(base, s, step, states, q0, accept, k0, J), dps - 15), flush=True)
    if not only or "bd1.9" in only:
        S = mpf(0)
        for n in range(1, 2 ** 9 + 1):
            if "11" not in bin(n):
                S += mpf(2) ** -n
        print("bd1.9", nstr(S, dps - 15))
