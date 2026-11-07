"""Experiment 192 -- SIAM 100-Digit Challenge, App. D Problem 17 (Python column).

Third Dirichlet eigenfunction of the Laplacian on the L-shaped domain and the
length of its zero set.  scipy/mpmath have no PDE eigensolver on a region
(a hand-built finite-difference Laplacian + eigsh would be a hand-rolled
discretisation, not a built-in), so every row is a SKIP.
See lshape_nodal_line.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

skip("d17.machine", "NDEigensystem")
skip("d17.machine-reduced", "NDEigenvalues")
skip("d17.wp50", "NDEigensystem")
