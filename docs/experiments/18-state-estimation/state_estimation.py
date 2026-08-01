#!/usr/bin/env python3
"""Experiment 18 -- State estimation: small dense matrices in a long loop.

Runs the same algorithm as ``state_estimation.m``, in the same order, with the
same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 state_estimation.py

WHAT IT MEASURES.  Every other benchmark in the suite is dominated by the cost
of ARRAYS.  This one is dominated by the cost of a CALL: 20000 iterations of
about ten operations on 6x6 and 2x2 matrices, for a total flop count a single
``dgemm`` would do in one call.

NumPy is not fast here either -- it is a Python loop around small BLAS calls,
which is exactly the regime the row is meant to expose.  It is nevertheless
the fastest of the three, and the ``README.md`` roadmap explains why and what
it would take for Mathilda to lead.

DETERMINISM.  The measurement sequence is a fixed pair of sinusoids, so
``trace(P)`` at the end is an exact cross-system check.  The filter is stable,
so P converges to a steady state rather than drifting.
"""

import math
import time

import numpy as np

# ---- shared reporting helpers (identical in every experiment file) --------


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs."""
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


# ---- model ---------------------------------------------------------------

KFN = 20000                # measurements
KFDT = 0.01                # sample interval

# Constant-acceleration transition in 2D: state is {x, y, vx, vy, ax, ay}.
kfF = np.zeros((6, 6))
for i in range(6):
    kfF[i, i] = 1.0
    if i < 4:
        kfF[i, i + 2] = KFDT
    if i < 2:
        kfF[i, i + 4] = 0.5 * KFDT ** 2
kfFT = kfF.T.copy()

kfH = np.zeros((2, 6))                 # observe x, y
kfH[0, 0] = 1.0
kfH[1, 1] = 1.0
kfHT = kfH.T.copy()
kfQ = 0.001 * np.eye(6)                # process noise
kfR = 0.05 * np.eye(2)                 # measurement noise
kfI = np.eye(6)

kfx0 = np.zeros(6)
kfP0 = np.eye(6)


def kfz(k):
    """A deterministic 'trajectory plus sensor noise'."""
    return np.array([math.sin(0.01 * k) + 0.02 * math.sin(7.1 * k),
                     math.cos(0.013 * k) + 0.02 * math.cos(5.3 * k)])


# ---- 1. the scalar filter ------------------------------------------------


def kfstep(st, k):
    """Textbook predict/update.  Every matrix is 6x6 or smaller."""
    xp = kfF @ st[0]                       # predict state
    pp = kfF @ st[1] @ kfFT + kfQ          # predict covariance
    yy = kfz(k) - kfH @ xp                 # innovation
    ss = kfH @ pp @ kfHT + kfR             # innovation covariance, 2x2
    kk = pp @ kfHT @ np.linalg.inv(ss)     # Kalman gain, 6x2
    return [xp + kk @ yy, (kfI - kk @ kfH) @ pp]


def kalman(m):
    st = [kfx0, kfP0]
    k = 1
    while k <= m:
        st = kfstep(st, k)
        k += 1
    return float(np.trace(st[1]))


# ---- 2. the ensemble filter ----------------------------------------------

ENN = 4096                 # ensemble members
ENS = 200                  # steps

enE0 = 0.01 * np.sin(np.arange(1, ENN + 1)[:, None] + 6 * np.arange(1, 7)[None, :])


def enstep(e, k):
    """The same estimator, array-shaped: the covariance is estimated from the
    spread of the ensemble rather than propagated from a model.  One 4096x6
    GEMM replaces the 6x6 algebra, so this row lives in the opposite regime."""
    pr = e @ kfFT                          # advance every member
    mn = pr.sum(0) / ENN                   # ensemble mean
    an = pr - mn[None, :]                  # anomalies
    pc = an.T @ an / (ENN - 1.0)           # sample covariance, 6x6
    kk = pc @ kfHT @ np.linalg.inv(kfH @ pc @ kfHT + kfR)
    inn = np.tile(kfz(k), (ENN, 1)) - pr @ kfHT
    return pr + inn @ kk.T


def enkf(m):
    e = enE0
    k = 1
    while k <= m:
        e = enstep(e, k)
        k += 1
    return float((e.sum(0) / ENN).sum())


def main():
    print("Experiment 18 -- state estimation")
    print("")

    bench("Kalman filter, 6 states, 20000 steps", lambda: kalman(KFN))
    check("Kalman filter (2000 steps)", kalman(2000))

    bench("ensemble Kalman, 4096 members, 200 steps", lambda: enkf(ENS))
    check("ensemble Kalman (20 steps)", enkf(20))

    # What one small dense operation costs.  NumPy pays Python-call overhead
    # here rather than kernel time, which is the honest comparison for the
    # .m file's numbers: both systems are measuring dispatch, not flops.
    print("")
    print("-- what one small dense operation costs --")
    sm6 = np.full((6, 6), 0.3) + 1.7 * np.eye(6)
    sv6 = np.arange(1.0, 7.0)
    sm2 = np.array([[2.0, 0.3], [0.3, 1.7]])
    sm60 = np.full((60, 60), 0.3) + 59.7 * np.eye(60)

    def rep(f, n):
        return lambda: [f() for _ in range(n)]

    bench("inv(2x2) x 20000", rep(lambda: np.linalg.inv(sm2), 20000))
    bench("inv(6x6) x 2000", rep(lambda: np.linalg.inv(sm6), 2000))
    bench("m6 @ v6  x 50000", rep(lambda: sm6 @ sv6, 50000))
    bench("m6 @ m6  x 20000", rep(lambda: sm6 @ sm6, 20000))
    bench("inv(60x60) x 2000", rep(lambda: np.linalg.inv(sm60), 2000))


if __name__ == "__main__":
    main()
