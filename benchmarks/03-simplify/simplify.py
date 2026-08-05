#!/usr/bin/env python3
"""Experiment 03 -- Simplify and FullSimplify (sympy column).

Same six kernels as ``simplify.m``, same order.

WHY SIMPLIFY IS HARD TO BENCHMARK HONESTLY.  Simplify has no unique right
answer, only a smaller one, so a timing alone rewards whichever system gives up
soonest.  Every row therefore checks a NUMERIC value at a rational point: any
correct simplification preserves it, and a system that returned its input
unchanged still passes the check while showing its cost.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy

from harness import bench, check, require

require(["sympy", "sympy:simplify", "sympy:cancel", "sympy:radsimp",
         "sympy:expand_trig"])

x = sympy.Symbol("x")
R = sympy.Rational


def at(e, v):
    return int(sympy.floor(sympy.re(sympy.N(e.subs(x, v) * 10 ** 6, 50))
                           + R(1, 2)))


t1 = sympy.sin(x) ** 2 + sympy.cos(x) ** 2
bench("Simplify trig Pythagorean", lambda: sympy.simplify(t1))
check("Simplify trig Pythagorean", at(sympy.simplify(t1), R(1, 3)))

t2 = 8 * sympy.sin(x) ** 4 - 8 * sympy.sin(x) ** 2 + 1
bench("Simplify quartic-to-Cos[4x]", lambda: sympy.simplify(t2))
check("Simplify quartic-to-Cos[4x]", at(sympy.simplify(t2), R(1, 3)))

t3 = (x ** 3 - 1) / (x - 1)
bench("Cancel removable singularity", lambda: sympy.cancel(t3))
check("Cancel removable singularity", at(sympy.cancel(t3), R(5, 2)))

t4 = sympy.sqrt(3 + 2 * sympy.sqrt(2))
bench("FullSimplify nested radical", lambda: sympy.simplify(sympy.radsimp(t4)))
check("FullSimplify nested radical",
      int(sympy.floor(sympy.N(sympy.simplify(sympy.radsimp(t4)) * 10 ** 6, 50)
                      + R(1, 2))))

t5 = sympy.log(sympy.exp(x)) + sympy.exp(sympy.log(x))
bench("Simplify log-exp collapse", lambda: sympy.simplify(t5))
check("Simplify log-exp collapse", at(sympy.simplify(t5), R(7, 4)))

t6 = sympy.sin(x) * sympy.sin(2 * x) * sympy.sin(3 * x) * sympy.sin(4 * x)
# TrigReduce's job is linearisation: products of sines -> a sum of cosines.
bench("TrigReduce product of 4 sines",
      lambda: sympy.simplify(sympy.expand(t6, trig=True)))
check("TrigReduce product of 4 sines",
      at(sympy.simplify(sympy.expand(t6, trig=True)), R(1, 3)))
