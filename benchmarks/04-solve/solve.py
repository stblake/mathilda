#!/usr/bin/env python3
"""Experiment 04 -- Solve: polynomial, systems, transcendental (sympy column).

Same six kernels as ``solve.m``, same order.

VALUE CHECKS ARE COUNTS AND SUMS, not solution sets: two systems legitimately
order roots differently and spell radicals differently, so the check is the
number of solutions plus the rounded sum of the roots -- both invariant under
ordering and under algebraic form.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import dense_upoly

require(["sympy", "sympy:solve", "sympy:nroots", "sympy:solveset"])

x, y, z = sympy.symbols("x y z")

e1 = sympy.prod([x - i for i in range(1, 7)])
bench("Solve factorable sextic", lambda: sympy.solve(sympy.Eq(e1, 0), x))
check("Solve factorable sextic", len(sympy.solve(sympy.Eq(e1, 0), x)))

e2 = x ** 4 - 10 * x ** 2 + 1
bench("Solve quartic with radicals", lambda: sympy.solve(sympy.Eq(e2, 0), x))
check("Solve quartic with radicals", len(sympy.solve(sympy.Eq(e2, 0), x)))

sys2 = [sympy.Eq(x ** 2 + y ** 2, 25), sympy.Eq(x - y, 1)]
bench("Solve 2x2 polynomial system", lambda: sympy.solve(sys2, [x, y]))
check("Solve 2x2 polynomial system", len(sympy.solve(sys2, [x, y])))

sys3 = [sympy.Eq(x + y + z, 6), sympy.Eq(x - y, 1), sympy.Eq(y - z, 1)]
# sympy returns a single dict for a determined linear system; Wolfram returns a
# list of one rule-list. Both have "one solution", so the check is 1 either way.
bench("Solve 3x3 linear system", lambda: sympy.solve(sys3, [x, y, z]))
_s3 = sympy.solve(sys3, [x, y, z])
check("Solve 3x3 linear system", 1 if isinstance(_s3, dict) else len(_s3))

p40 = sympy.Poly(dense_upoly(40, x), x)
bench("NSolve degree 40", lambda: p40.nroots(n=15, maxsteps=200))
check("NSolve degree 40", len(p40.nroots(n=15, maxsteps=200)))

bench("Solve sextic, root sum", lambda: sum(sympy.solve(sympy.Eq(e1, 0), x)))
check("Solve sextic, root sum", int(sum(sympy.solve(sympy.Eq(e1, 0), x))))
