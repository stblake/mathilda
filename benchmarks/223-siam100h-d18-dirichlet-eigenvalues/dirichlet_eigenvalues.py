"""Experiment 223 -- SIAM 100-Digit Challenge holdout A, family hd18 (Python column).

Dirichlet eigenvalues of -Laplacian on planar regions (see
dirichlet_eigenvalues.m).  scipy and mpmath have no eigensolver for a PDE on
a region (no FEM): every row is a SKIP, as in exp 193."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

for k in range(1, 11):
    for tier in ("machine", "wp50"):
        skip("hd18.%d %s" % (k, tier), "NDEigenvalues (FEM)")
