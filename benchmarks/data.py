#!/usr/bin/env python3
"""benchmarks/data.py -- shared, deterministic input generators.

The Python twin of ``data.m``: every function here builds the **same
mathematical object** its ``.m`` counterpart does, so a check value is
comparable across the three systems.

WHY THIS IS SHARED.  The value gate is only meaningful if both halves of an
experiment compute over the same inputs.  One copy per language makes that true
by construction rather than by two authors agreeing on a formula.

Everything is exactly reproducible -- integer coefficients, closed forms, no
floating-point randomness.  Where an experiment needs random float data for its
*timing*, it draws that locally and checks a separate small deterministic case:
the three systems cannot be made to draw the same float stream, so a float sum
over random data is never a check.
"""

# sympy is imported lazily so the numeric experiments (groups B, C, D) do not
# pay sympy's ~1 s import cost just to reach a numpy kernel.


# ---- univariate polynomials ---------------------------------------------


def dense_upoly(n, x):
    """Dense degree-n polynomial, small reproducible integer coefficients.

    Matches ``denseUPoly`` in data.m: sum of ((k^2+3k+1) mod 17) x^k, plus 1.
    Not random -- mild coefficient growth keeps Factor's cost about degree
    rather than about coefficient size.
    """
    return sum(((k * k + 3 * k + 1) % 17) * x ** k for k in range(n + 1)) + 1


def factorable_upoly(k, x):
    """Product of distinct quadratics -- a factoring target with a known answer."""
    out = 1
    for i in range(1, k + 1):
        out = out * (x ** 2 + i * x + (i + 1))
    return out


def sparse_upoly(n, x):
    """Sparse high degree: cheap to store, expensive to factor."""
    return x ** n + x ** (n // 2) + 1


# ---- rational functions -------------------------------------------------


def rat_fn(n, x):
    """Proper rational function, squarefree degree-n denominator."""
    den = 1
    for i in range(1, n + 1):
        den = den * (x - i)
    return dense_upoly(n - 1, x) / den


def rat_fn_repeated(n, x):
    """Repeated denominator factors -- the branch Hermite reduction is for."""
    den = 1
    for i in range(1, n + 1):
        den = den * (x - i) ** 2
    return dense_upoly(n, x) / den


# ---- multivariate -------------------------------------------------------


def multi_poly(k, variables):
    import sympy
    s = sum(variables)
    out = 1
    for i in range(1, k + 1):
        out = out * (s + i)
    return sympy.expand(out)


# ---- polynomial systems (Groebner targets) ------------------------------


def cyclic_system(variables):
    """cyclic-n, the standard Groebner stress family. Matches cyclicSystem."""
    n = len(variables)
    eqs = []
    for d in range(1, n):
        total = 0
        for i in range(n):
            term = 1
            for k in range(d):
                term = term * variables[(i + k) % n]
            total = total + term
        eqs.append(total)
    prod = 1
    for v in variables:
        prod = prod * v
    eqs.append(prod - 1)
    return eqs


# ---- matrices -----------------------------------------------------------


def spd_int_matrix(n):
    """Symmetric positive-definite, exactly reproducible. Matches spdIntMatrix.

    Returned as float64 for the numeric column; the ``.m`` side holds it
    exactly, so any check must be a rounded scalar rather than a bit pattern.
    """
    import numpy as np
    m = np.empty((n, n), dtype=np.float64)
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            m[i - 1, j - 1] = (n + 2) if i == j else 1.0 / (1 + abs(i - j))
    return m


def hilbert_like(n):
    import numpy as np
    i = np.arange(1, n + 1)
    return 1.0 / (i[:, None] + i[None, :] - 1)


# ---- integer data -------------------------------------------------------


def lcg_ints(n, modulus):
    """Reproducible integers independent of either system's RNG stream.

    A linear congruential sequence, identical to ``lcgInts`` in data.m.  This is
    what makes a Sort or a Total check comparable across three systems -- the
    built-in generators cannot be aligned, so we bring our own.
    """
    out = []
    s = 12345
    for _ in range(n):
        s = (1103515245 * s + 12345) % 2147483648
        out.append(s % modulus + 1)
    return out
