"""Experiment 252 -- SIAM 100-digit SEALED HOLDOUT B, family D17 (Python column).

The nodal set of the k-th Dirichlet eigenfunction on a region (see
nodal_lines.m).  scipy and mpmath have no FEM eigensolver on a geometric
region: every row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip, require

require(["scipy.sparse.linalg:eigsh"])
for k in range(1, 11):
    skip("bd17.%d machine" % k, "NDEigensystem")
    skip("bd17.%d wp50" % k, "NDEigensystem")
