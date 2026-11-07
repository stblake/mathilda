"""Experiment 189 -- SIAM 100-Digit Challenge, App. D Problem 14 (Python column).

The idiomatic statement is a finite-element solve of a mixed boundary-value
problem on a pentagon (see pentagon_poisson.m).  Neither scipy nor mpmath has
a PDE / finite-element solver, so both rows are SKIP rather than a hand-written
discretisation."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

skip("d14.machine", "NDSolveValue (FEM)")
skip("d14.wp50", "NDSolveValue (FEM)")
