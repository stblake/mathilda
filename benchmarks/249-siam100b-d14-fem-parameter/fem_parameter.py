"""Experiment 249 -- SIAM 100-Digit Challenge holdout B for App. D Problem 14 (Python column).

Every case is a finite-element boundary-value solve inside a root finder (see
fem_parameter.m).  Neither scipy nor mpmath has a PDE / finite-element solver,
so all rows are SKIP rather than a hand-written discretisation (as in exp. 189)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

for tag in ("machine", "wp50"):
    for k in range(1, 11):
        skip("bd14.%d %s" % (k, tag), "NDSolveValue (FEM)")
