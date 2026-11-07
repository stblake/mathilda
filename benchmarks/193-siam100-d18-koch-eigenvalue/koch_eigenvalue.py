"""Experiment 193 -- SIAM 100-Digit Challenge, App. D Problem 18 (Python column).

lambda_1 of the Dirichlet Laplacian on the Koch snowflake.  scipy/mpmath have
no PDE eigensolver on a planar region (a hand-built FEM/FD matrix would be a
hand-rolled discretisation, not a built-in): every row is a SKIP.
See koch_eigenvalue.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip

skip("d18.machine", "NDEigenvalues")
skip("d18.machine level=0", "NDEigenvalues")
skip("d18.wp50", "NDEigenvalues")
