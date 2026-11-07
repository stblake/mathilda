"""Experiment 221 -- SIAM 100-Digit Challenge holdout A, family hd16 (Python column).

Heat-equation crossing times (see heat_crossings.m).  The problem is an
initial-boundary-value problem; scipy has no PDE solver and mpmath has none,
so every row is a SKIP (as for the PDE rows of exps 173 and 191)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

for k in range(1, 11):
    for tier in ("machine", "wp50"):
        skip("hd16.%d %s" % (k, tier), "PDE solver")
