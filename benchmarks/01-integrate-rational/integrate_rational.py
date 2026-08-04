#!/usr/bin/env python3
"""Experiment 01 -- Symbolic integration of rational functions (sympy column).

The same six kernels as ``integrate_rational.m``, in the same order.

    python3 integrate_rational.py

WHAT IT MEASURES.  sympy's ``integrate`` on the rational-function branch:
partial fractions, Hermite reduction for repeated factors, and the logarithmic
part.  The Mathilda side is ``src/calculus/intrat.c`` and ``src/parfrac.c``.

HOW THE VALUE CHECK SURVIVES A DIFFERENT ANSWER.  Two systems can return
algebraically different antiderivatives that are both correct -- they differ by
a constant and by how the logarithmic part is grouped.  So the check is the
DEFINITE integral recovered from the antiderivative, ``F(b) - F(a)``, which is
invariant under both.  The interval sits clear of every pole, so no branch cut
is crossed.

NOT A NUMPY ROW.  sympy is pure Python; this column is "what a Python CAS
does", not "what C can do".  That is the right comparison here -- Mathilda is
also a CAS doing symbolic work -- but it means a Mathilda win on these rows is
a weaker claim than a win on a numeric row against a BLAS call.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import dense_upoly, rat_fn, rat_fn_repeated

require(["sympy", "sympy:integrate", "sympy:apart", "sympy:together"])

x = sympy.Symbol("x")


def defint(expr, a, b):
    """F(b) - F(a) from an antiderivative, rounded to an integer at 1e6.

    Matches the .m side's Round[10^6 (F[b] - F[a])].  Uses exact rational
    endpoints so the two systems substitute the same number, then evaluates to
    50 digits before rounding -- rounding a float64 at 1e6 would make the check
    sensitive to the last bit.
    """
    F = sympy.integrate(expr, x)
    d = F.subs(x, sympy.Rational(b)) - F.subs(x, sympy.Rational(a))
    return int(sympy.floor(sympy.re(sympy.N(d * 10 ** 6, 50)) + sympy.Rational(1, 2)))


# ---- 1. Squarefree denominator: partial fractions + log part -------------

f6 = rat_fn(6, x)
bench("integrate rational, squarefree deg 6", lambda: sympy.integrate(f6, x))
check("integrate rational, squarefree deg 6",
      defint(f6, sympy.Rational(13, 2), sympy.Rational(15, 2)))

# ---- 2. Repeated denominator factors: the Hermite reduction branch -------

f4r = rat_fn_repeated(4, x)
bench("integrate rational, repeated deg 4", lambda: sympy.integrate(f4r, x))
check("integrate rational, repeated deg 4",
      defint(f4r, sympy.Rational(9, 2), sympy.Rational(11, 2)))

# ---- 3. Irreducible denominator: no rational part, all logarithmic -------

# x**4 + 1 factors over Q into two irreducible quadratics, so the log part needs
# a real algebraic extension, and the antiderivative is closed form in log/atan.
#
# NOT x**5 + x + 1, the first choice: its antiderivative carries an unresolved
# root sum that the two systems render differently, which broke the value check
# (not the timing). A row that cannot be checked cannot be trusted.
fq = 1 / (x ** 4 + 1)
bench("integrate rational, irreducible quartic", lambda: sympy.integrate(fq, x))
check("integrate rational, irreducible quartic",
      defint(fq, sympy.Rational(3, 2), sympy.Rational(5, 2)))

# ---- 4. Polynomial integration: the trivial baseline --------------------

p40 = dense_upoly(40, x)
bench("integrate polynomial, deg 40", lambda: sympy.integrate(p40, x))
check("integrate polynomial, deg 40",
      int(sympy.floor(sympy.N(sympy.integrate(p40, x).subs(x, 1) * 10 ** 6, 50)
                      + sympy.Rational(1, 2))))

# ---- 5. Partial fractions alone, without the integration ----------------

f8 = rat_fn(8, x)
bench("Apart, squarefree deg 8", lambda: sympy.apart(f8, x))
check("Apart, squarefree deg 8",
      int(sympy.floor(sympy.N(sympy.apart(f8, x).subs(x, sympy.Rational(19, 2))
                              * 10 ** 6, 50) + sympy.Rational(1, 2))))

# ---- 6. Together: the inverse direction ---------------------------------

summands = sum(1 / (x - i) for i in range(1, 9))
bench("Together, 8 summands", lambda: sympy.together(summands))
check("Together, 8 summands",
      int(sympy.floor(sympy.N(sympy.together(summands).subs(
          x, sympy.Rational(19, 2)) * 10 ** 6, 50) + sympy.Rational(1, 2))))
