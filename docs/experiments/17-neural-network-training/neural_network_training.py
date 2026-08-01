#!/usr/bin/env python3
"""Experiment 17 -- Neural-network training: a chain of transposed GEMMs.

Runs the same algorithm as ``neural_network_training.m``, in the same order,
with the same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 neural_network_training.py

WHAT IT MEASURES.  The third sweep's logistic-regression row has ONE weight
matrix.  Real training has a chain, and the backward pass is a chain of
TRANSPOSED products -- ``a1.T @ d2``, ``d2 @ W2.T``, ``X.T @ d1``.  A
two-layer MLP is the smallest kernel that exercises the whole chain.

WHERE NUMPY HAS A STRUCTURAL ADVANTAGE, STATED PLAINLY.  ``X.T`` is a *view*:
it costs nothing and BLAS is then told to read the buffer transposed.  Both
CAS copy.  ``X.T`` in this file is therefore free and the corresponding
``Transpose[nnX]`` in the ``.m`` file is 2 ms a step, recomputed 100 times
because it is loop-invariant but re-evaluated.  That is a capability
difference, not a tuning difference, and it is the largest single item on this
experiment's roadmap.

BIASES ARE A COLUMN OF X, NOT A SEPARATE VECTOR -- what a performance-minded
implementation does in all three languages, and it keeps the three columns
running identical arithmetic.

DETERMINISM.  The timed run uses random data.  The check trains a
deterministic 256 x 17 -> 8 -> 4 network for 20 steps and compares the
cross-entropy loss, which depends on every weight and gradient along the way.
"""

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


def nest(f, x, m):
    for _ in range(m):
        x = f(x)
    return x


def relu(z):
    return np.maximum(z, 0.0)


# ---- the timed network ---------------------------------------------------

NNB = 1024                 # minibatch
NNI = 785                  # inputs, including the constant bias column
NNH = 128                  # hidden units
NNO = 10                   # classes
NNLR = 0.05                # learning rate

nnX = np.concatenate([np.random.uniform(-1, 1, (NNB, NNI - 1)),
                      np.ones((NNB, 1))], 1)
nnY = np.zeros((NNB, NNO))
nnY[np.arange(NNB), (np.arange(1, NNB + 1) % NNO)] = 1.0

# Deterministic initial weights: no RNG, so the two CAS start identically.
nnW1 = 0.05 * np.sin(np.arange(1, NNI + 1)[:, None] + NNI * np.arange(1, NNH + 1)[None, :])
nnW2 = 0.05 * np.cos(np.arange(1, NNH + 1)[:, None] + NNH * np.arange(1, NNO + 1)[None, :])


def sfmax(z):
    """Row-wise softmax, shifted by the row maximum for numerical stability."""
    ex = np.exp(z - z.max(1)[:, None])
    return ex / ex.sum(1)[:, None]


def mlpstep(w):
    """One SGD step: forward, softmax cross-entropy, backward, update."""
    z1 = nnX @ w[0]                        # 1024x785 @ 785x128
    a1 = relu(z1)
    z2 = a1 @ w[1]                         # 1024x128 @ 128x10
    sm = sfmax(z2)
    d2 = (sm - nnY) / NNB                  # dL/dz2 for cross-entropy
    d1 = (d2 @ w[1].T) * (z1 >= 0)         # backprop through ReLU
    g2 = a1.T @ d2
    g1 = nnX.T @ d1                        # nnX.T is a VIEW here, free
    return [w[0] - NNLR * g1, w[1] - NNLR * g2]


def mlptrain(m):
    return nest(mlpstep, [nnW1, nnW2], m)


# ---- the deterministic check network -------------------------------------

CB, CI, CH, CO, CLR = 256, 17, 8, 4, 0.1
cX = np.concatenate([np.sin(0.3 * (np.arange(1, CB + 1)[:, None]
                                   + CI * np.arange(1, CI)[None, :])),
                     np.ones((CB, 1))], 1)
cY = np.zeros((CB, CO))
cY[np.arange(CB), (np.arange(1, CB + 1) % CO)] = 1.0
cW1 = 0.2 * np.sin(np.arange(1, CI + 1)[:, None] + CI * np.arange(1, CH + 1)[None, :])
cW2 = 0.2 * np.cos(np.arange(1, CH + 1)[:, None] + CH * np.arange(1, CO + 1)[None, :])


def cstep(w):
    z1 = cX @ w[0]
    a1 = relu(z1)
    z2 = a1 @ w[1]
    sm = sfmax(z2)
    d2 = (sm - cY) / CB
    d1 = (d2 @ w[1].T) * (z1 >= 0)
    g2 = a1.T @ d2
    g1 = cX.T @ d1
    return [w[0] - CLR * g1, w[1] - CLR * g2]


def closs(m):
    w = nest(cstep, [cW1, cW2], m)
    sm = sfmax(relu(cX @ w[0]) @ w[1])
    return float(-(cY * np.log(sm + 1e-12)).sum() / CB)


# ---- inference -----------------------------------------------------------

NNQ = 8192
nnXq = np.concatenate([np.random.uniform(-1, 1, (NNQ, NNI - 1)),
                       np.ones((NNQ, 1))], 1)


def mlpinfer():
    return sfmax(relu(nnXq @ nnW1) @ nnW2).sum(1)


def main():
    print("Experiment 17 -- neural-network training")
    print("")

    bench("MLP training, 785-128-10, batch 1024, 100 steps",
          lambda: mlptrain(100))
    check("MLP training (loss after 20 deterministic steps)", closs(20))

    bench("MLP inference, 8192 x 785 forward pass", mlpinfer)
    check("MLP inference (loss after 1 step)", closs(1))

    # One training step, split.  The two lines that dominate the .m file --
    # the loop-invariant transpose and the mixed-dtype mask -- are free and
    # nearly free here, which is what makes this the useful control.
    print("")
    print("-- one training step, split (per step) --")
    sz1 = nnX @ nnW1
    sa1 = relu(sz1)
    sd2 = np.random.rand(NNB, NNO)
    sda = sd2 @ nnW2.T
    sdz = sda * (sz1 >= 0)

    bench("nnX @ W1        (forward GEMM)", lambda: nnX @ nnW1)
    bench("maximum(z1, 0)  (ReLU)", lambda: relu(sz1))
    bench("z2.max(1)       (softmax shift)", lambda: (sa1 @ nnW2).max(1))
    bench("nnX.T           (a VIEW, free)", lambda: nnX.T)
    bench("a1.T", lambda: sa1.T)
    bench("a1.T @ d2", lambda: sa1.T @ sd2)
    bench("da1 * (z1 >= 0) (mask)", lambda: sda * (sz1 >= 0))
    bench("nnX.T @ dz1", lambda: nnX.T @ sdz)
    bench("W1 - lr g1      (update)", lambda: nnW1 - NNLR * (nnX.T @ sdz))


if __name__ == "__main__":
    main()
