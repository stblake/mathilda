"""Experiment 191 -- SIAM 100-Digit Challenge, App. D Problem 16 (Python column).

The problem starts from the Daubechies scaling function phi_2.  scipy 1.17
has no Daubechies scaling function (scipy.signal.daub / cascade were removed
in scipy 1.15) and mpmath has none; PyWavelets is not part of the reference
stack.  Every row is a SKIP.  See daubechies_heat.m."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siam100_common import skip, require

require(["scipy.signal:cascade", "mpmath:quad"])

skip("d16.machine", "DaubechiesScalingFunction")
skip("d16.machine-reduced", "DaubechiesScalingFunction")
skip("d16.wp50", "DaubechiesScalingFunction")
