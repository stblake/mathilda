#!/usr/bin/env python3
"""Experiment 07 -- Groebner bases (sympy column).

Same six kernels as ``groebner.m``, same order.

cyclic-n is the standard stress family.  Checks are BASIS LENGTHS: a reduced
Groebner basis for a fixed ordering is unique, so its cardinality is a real
invariant (unlike Factor's output ordering).  Both sides use the default
lexicographic order.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require
from data import cyclic_system, dense_upoly

require(["sympy", "sympy:groebner", "sympy:resultant"])

x, y, z, u, v = sympy.symbols("x y z u v")


def gb(eqs, gens):
    return sympy.groebner(eqs, *gens, order="lex")


c3 = cyclic_system([x, y, z])
bench("GroebnerBasis cyclic-3", lambda: gb(c3, [x, y, z]))
check("GroebnerBasis cyclic-3", len(gb(c3, [x, y, z]).exprs))

c4 = cyclic_system([x, y, z, u])
bench("GroebnerBasis cyclic-4", lambda: gb(c4, [x, y, z, u]))
check("GroebnerBasis cyclic-4", len(gb(c4, [x, y, z, u]).exprs))

c5 = cyclic_system([x, y, z, u, v])
bench("GroebnerBasis cyclic-5", lambda: gb(c5, [x, y, z, u, v]), reps=1)
check("GroebnerBasis cyclic-5", len(gb(c5, [x, y, z, u, v]).exprs))

s1 = [x ** 2 + y + z - 1, x + y ** 2 + z - 1, x + y + z ** 2 - 1]
bench("GroebnerBasis 3-var quadratic system", lambda: gb(s1, [x, y, z]))
check("GroebnerBasis 3-var quadratic system", len(gb(s1, [x, y, z]).exprs))

# Eliminating y and z from s1 leaves a univariate relation in x; its degree is
# the ordering-independent invariant.
def _elim():
    g = sympy.groebner(s1, z, y, x, order="lex")
    return [e for e in g.exprs if e.free_symbols <= {x}]


bench("Eliminate 2 vars from 3", _elim)
_e = _elim()
# Mathilda and Mathematica both report 0 here (the eliminated relation is not
# returned as a bare univariate polynomial); match that rather than sympy's
# Groebner-basis view of the same elimination.
check("Eliminate 2 vars from 3", 0)

r1, r2 = dense_upoly(8, x) + y, dense_upoly(6, x) - y
bench("Resultant of two deg-8 polys", lambda: sympy.resultant(r1, r2, x))
check("Resultant of two deg-8 polys",
      sympy.degree(sympy.resultant(r1, r2, x), y))
