#!/usr/bin/env python3
"""Experiment 15 -- Option pricing: trees, PDEs, and the positive part.

Runs the same algorithm as ``option_pricing.m``, in the same order, with the
same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 option_pricing.py

WHAT IT MEASURES.  Two things no other benchmark does: a working vector that
SHRINKS (a binomial tree goes from 4001 elements to 1, crossing the packing
threshold inside the loop), and an early-exercise projection --
``max(continuation, payoff)``, which is the positive part and the same
operation as a rectified linear unit.

TWO PRICERS, ON PURPOSE.  A 4000-step Cox-Ross-Rubinstein tree and a
1000 x 25000 explicit finite-difference solver in log-price.  They are
independent discretisations of the same problem and agree to about 0.6%, which
is a stronger statement than either agreeing with itself across three systems.

THE FD STABILITY CONDITION IS LOAD-BEARING.  An explicit scheme needs
``sigma^2 dt / dx^2 <= 1/2``.  At 1000 x 25000 that ratio is 0.433.  The first
draft used 2000 x 7000, where it is 6.19; the scheme diverged, every system
returned NaN or a plausible-looking wrong number, and only the cross-system
value check caught it.
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


def relu(z):
    """The positive part.  ``np.maximum(z, 0)`` is the idiomatic NumPy
    spelling; the ``.m`` file uses ``Ramp``, which this experiment caused to
    exist in Mathilda."""
    return np.maximum(z, 0.0)


# ---- contract ------------------------------------------------------------

BTN = 4000                 # tree steps
BTT = 1.0                  # years to expiry
BTR = 0.03                 # risk-free rate
BTSIG = 0.25               # volatility
BTK = 100.0                # strike
BTS0 = 100.0               # spot

# ---- 1. Cox-Ross-Rubinstein binomial tree, American put ------------------

BTDT = BTT / BTN
BTU = math.exp(BTSIG * math.sqrt(BTDT))                       # up factor
BTP = (math.exp(BTR * BTDT) - 1.0 / BTU) / (BTU - 1.0 / BTU)  # risk-neutral
BTQ = 1.0 - BTP
BTDC = math.exp(-BTR * BTDT)                                  # discount


def amtree():
    """Backward induction.

    At level k the spot lattice is ``S0 u^(2j-k)``; stepping back multiplies
    the surviving nodes by ``u``, so the spot vector is maintained by
    ``BTU * s[:-1]`` rather than rebuilt.  The vector shrinks by one element
    per iteration, which is the point of this row.
    """
    s = BTS0 * BTU ** (2.0 * np.arange(0, BTN + 1) - BTN)
    v = relu(BTK - s)                              # terminal payoff
    k = BTN
    while k > 0:
        s = BTU * s[:-1]
        v = BTDC * (BTP * v[1:] + BTQ * v[:-1])    # continuation value
        v = np.maximum(v, relu(BTK - s))           # early exercise
        k -= 1
    return float(v[0])


# ---- 2. explicit finite difference in log-price, American put ------------

FDM = 1000                 # space points
FDN = 25000                # time steps
FDXL = 1.2                 # log-price half-width

FDDX = 2.0 * FDXL / (FDM - 1)
FDDT = 1.0 / FDN
fdx = -FDXL + np.arange(0, FDM) * FDDX
fdpay = relu(BTK - BTS0 * np.exp(fdx))             # payoff, once

# Constant coefficients: the log transform removes the S-dependence, so the
# three stencil weights are scalars and the sweep is one shifted triple.
FDA = 0.5 * BTSIG ** 2 * FDDT / FDDX ** 2 - 0.5 * (BTR - 0.5 * BTSIG ** 2) * FDDT / FDDX
FDC = 0.5 * BTSIG ** 2 * FDDT / FDDX ** 2 + 0.5 * (BTR - 0.5 * BTSIG ** 2) * FDDT / FDDX
FDB = 1.0 - BTSIG ** 2 * FDDT / FDDX ** 2 - BTR * FDDT

FDLO = BTK - BTS0 * math.exp(-FDXL)                # deep in the money
FDMID = (FDM + 1) // 2                             # the at-the-money index


def amfd():
    """Dirichlet boundaries by construction: the stencil is applied to the
    interior only and the two edges are re-attached with a concatenate."""
    v = fdpay.copy()
    k = 0
    while k < FDN:
        vi = FDA * v[:-2] + FDB * v[1:-1] + FDC * v[2:]
        v = np.concatenate(([FDLO], vi, [0.0]))
        v = np.maximum(v, fdpay)                   # early exercise
        k += 1
    return float(v[FDMID - 1])


# ---- 3. Monte-Carlo value at risk ----------------------------------------

VRN = 250000               # scenarios
VRK = 64                   # assets
vrs = np.random.uniform(-0.02, 0.02, (VRN, VRK))   # timed: random
vrw = np.full(VRK, 1.0 / VRK)                      # equal weights
vrcs = ((7 * np.arange(1, 20001)[:, None] + 3 * np.arange(1, VRK + 1)[None, :])
        % 41 - 20) / 1000.0


def mcvar(s, w, q):
    """Portfolio P&L is one dgemv; the risk number is a tail mean of the
    sorted result -- an order statistic, so this row is a Dot plus a Sort."""
    pnl = np.sort(s @ w)
    return float(pnl[:q].sum() / q)


def main():
    print("Experiment 15 -- option pricing")
    print("")

    bench("binomial American put, 4000 steps", amtree)
    check("binomial American put", amtree())

    bench("explicit FD American put, 1000 x 25000", amfd)
    check("explicit FD American put", amfd())

    bench("Monte-Carlo VaR, 250000 x 64", lambda: mcvar(vrs, vrw, 10000))
    check("Monte-Carlo VaR", mcvar(vrcs, vrw, 200))

    print("")
    print("the two pricers agree to about 0.6%, which is what says they are")
    print("pricing the same option rather than each reproducing itself")

    # The .m file compares three spellings of the positive part, because until
    # this experiment only one of them worked in Mathilda.  NumPy has only the
    # one idiomatic spelling, timed here as the reference.
    print("")
    print("-- the positive part, 10^6 float64 --")
    pv = np.random.uniform(-1, 1, 10 ** 6)
    bench("maximum(v, 0.)", lambda: np.maximum(pv, 0.0))
    bench("v * (v >= 0)", lambda: pv * (pv >= 0))


if __name__ == "__main__":
    main()
