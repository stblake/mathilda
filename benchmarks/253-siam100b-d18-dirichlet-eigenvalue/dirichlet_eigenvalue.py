"""Experiment 253 -- SIAM 100-digit SEALED HOLDOUT B, family D18 (Python column).

The smallest Dirichlet eigenvalue of -Laplacian on a region (see
dirichlet_eigenvalue.m).  scipy and mpmath have no PDE / FEM eigensolver on a
geometric region (scipy.sparse.linalg.eigsh needs a matrix the user would
have to discretise by hand): every row is a SKIP."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip, require

require(["scipy.sparse.linalg:eigsh"])
for k in range(1, 11):
    skip("bd18.%d machine" % k, "NDEigenvalues")
    skip("bd18.%d wp50" % k, "NDEigenvalues")
