"""benchmarks/siam100_common.py -- the Python twin of siam100_common.m.

Imported by every SIAM 100-digit experiment's .py half (scipy for the machine
rows, mpmath for the 50-digit rows)::

    import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from siam100_common import vrun, skip, require, mp

Emits the same lines as the .m side, with the same labels:

    BENCH  <label>  <ms>          only for an answer
    CHECK  <label>  {m, e}        sig(v, d), or timeout/unevaluated/failed/complex
    VALUE  <label>  {m, e}        every digit the result carries

sig(v, d) = {sign * leading d significant digits (rounded), decimal exponent},
computed in exact decimal arithmetic so a float prints its true binary value
and an mpf its full working precision.

A case is capped with SIGALRM: the alarm raises inside the integrand / the
objective the next time control returns to Python bytecode, which scipy and
mpmath both do on every function evaluation.
"""

import signal
import time
from decimal import Decimal, localcontext

import mpmath
from mpmath import mp

from harness import require, check   # noqa: F401  (re-exported)

CAP = 120


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def _digits_of(v):
    """Significant digits the value carries (17 for a double)."""
    if isinstance(v, mpmath.mpf):
        return max(1, int(mp.prec * 0.30103))
    return 17


def _to_decimal(v):
    if isinstance(v, mpmath.mpf):
        # nstr with enough digits is exact enough for scoring at <= prec.
        return Decimal(mpmath.nstr(v, _digits_of(v) + 5, strip_zeros=False,
                                   min_fixed=1, max_fixed=0))
    return Decimal(float(v))          # exact binary value of the double


def sig(v, d):
    """{sign * leading d significant digits, decimal exponent} as text."""
    x = _to_decimal(v)
    if x == 0:
        return "{0, 0}"
    with localcontext() as ctx:
        ctx.prec = d + 40
        a = abs(x)
        e = a.adjusted()                       # floor(log10 a)
        m = (a.scaleb(-(e - d + 1))).quantize(Decimal(1),
                                              rounding="ROUND_HALF_EVEN")
        if m >= Decimal(10) ** d:              # rounding carried a digit
            e += 1
            m = (a.scaleb(-(e - d + 1))).quantize(Decimal(1),
                                                  rounding="ROUND_HALF_EVEN")
        s = -1 if x < 0 else 1
    return "{%d, %d}" % (s * int(m), e)


def _real(v):
    """The real answer, or a string meaning 'not an answer'."""
    if v is None:
        return "unevaluated"
    if isinstance(v, str):
        return v
    if isinstance(v, (complex, mpmath.mpc)):
        if v.imag == 0:
            v = v.real
        else:
            return "complex"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "unevaluated"
    if f != f:
        return "failed"          # NaN
    if f in (float("inf"), float("-inf")) and not isinstance(v, mpmath.mpf):
        return "failed"
    return v


def vrun(label, fn, d, reps=1, cap=CAP, dps=None):
    """Time fn() (warm-up + min of reps when reps > 1) and emit the case.

    dps: the mpmath working precision for an MPFR row.  It is set around both
    the computation and the printing, because an mpf prints (and is scored) at
    the precision in force when it is formatted, not the one it was made at.
    """
    saved_dps = mp.dps
    if dps is not None:
        mp.dps = dps
    try:
        _vrun(label, fn, d, reps, cap)
    finally:
        mp.dps = saved_dps


def _vrun(label, fn, d, reps, cap):
    old = signal.signal(signal.SIGALRM, _alarm)
    ts, v = [], None
    try:
        signal.alarm(int(cap))
        try:
            if reps > 1:
                fn()
            for _ in range(reps):
                t0 = time.perf_counter()
                v = fn()
                ts.append(time.perf_counter() - t0)
        finally:
            signal.alarm(0)
    except _Timeout:
        v = "timeout"
    except Exception as e:                      # a solver raising = failed
        v = "failed"
        print("# %s: %s: %s" % (label, type(e).__name__, e), flush=True)
    finally:
        signal.signal(signal.SIGALRM, old)
    v = _real(v)
    if isinstance(v, str):
        check(label, v)
        return
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * min(ts)), flush=True)
    check(label, sig(v, d))
    print("VALUE\t%s\t%s" % (label, sig(v, _digits_of(v))), flush=True)


def skip(label, head):
    """A case this system cannot express with its built-in heads."""
    print("SKIP\t%s\t%s" % (label, head), flush=True)
