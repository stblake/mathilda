"""Experiment 222 -- SIAM 100-Digit Challenge holdout A, family hd17 (Python column).

Nodal-line lengths of Dirichlet eigenfunctions (see
nodal_lines.m).  scipy and mpmath have no eigensolver for a PDE on
a region (no FEM): every row is a SKIP, as in exp 192."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

for k in range(1, 11):
    for tier in ("machine", "wp50"):
        skip("hd17.%d %s" % (k, tier), "NDEigensystem (FEM)")
