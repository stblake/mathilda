#!/usr/bin/env python3
"""Experiment 13 -- Molecular dynamics: Lennard-Jones with a cut-off (NumPy).

Runs the same algorithm as ``molecular_dynamics.m``, in the same order, with
the same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 molecular_dynamics.py

WHAT IT MEASURES.  Lennard-Jones differs from plain N-body gravity in one
respect that matters to an array language: it multiplies through a
DATA-DEPENDENT mask.  The question the row answers is whether a comparison
result stays on the machine buffer or forces the 2048 x 2048 interaction
matrix back through the interpreter.

DETERMINISM.  The configuration is a deterministically perturbed simple-cubic
lattice, so the whole benchmark reproduces across systems rather than only its
check -- and the check is the potential energy, which every one of the 2
million pairs contributes to.

THE SELF-INTERACTION.  A large constant on the diagonal of ``r2`` makes
``i == j`` contribute nothing (1/r2 is 1e-12 rather than infinite, and the
cut-off mask is 0 there anyway).  That keeps the kernel branch-free, which is
the whole point of the array formulation.

WHY THERE IS NO SHARED HELPER.  ``ljf`` and ``lje`` recompute the same six
intermediates rather than sharing them.  Factoring them out is the obvious
thing to do and is a trap in Mathilda -- five of the six are float64 and the
mask is int64, and a mixed-dtype tuple of packed rows cannot be absorbed into
one buffer.  The ``.m`` file measures that; NumPy has no equivalent problem,
which is exactly what makes it a useful control.
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
    """Nest[f, x, m] -- apply f to x, m times."""
    for _ in range(m):
        x = f(x)
    return x


# ---- configuration -------------------------------------------------------

NMD = 2048                 # atoms
MDRC2 = 6.25               # squared cut-off radius, 2.5 sigma
MDDT = 0.0005              # Verlet time step

_i = np.arange(0, NMD, dtype=float)
mdx = (_i.astype(int) % 13).astype(float) + 0.13 * np.sin(1.7 * _i)
mdy = ((_i.astype(int) // 13) % 13).astype(float) + 0.13 * np.cos(2.3 * _i)
mdz = (_i.astype(int) // 169).astype(float) + 0.13 * np.sin(0.9 * _i)

mdbig = 1e12 * np.eye(NMD)             # kills the self-interaction
mdzero = np.zeros(NMD)


# ---- kernels -------------------------------------------------------------


def ljf(xs, ys, zs):
    """Lennard-Jones force on every atom.  Returns [fx, fy, fz]."""
    dx = xs[:, None] - xs[None, :]
    dy = ys[:, None] - ys[None, :]
    dz = zs[:, None] - zs[None, :]
    r2 = dx * dx + dy * dy + dz * dz + mdbig
    mk = (MDRC2 - r2 >= 0).astype(float)      # 1 inside the cut-off, else 0
    i2 = 1.0 / r2
    i6 = i2 * i2 * i2                         # (sigma/r)^6
    ff = (48.0 * i6 * i6 - 24.0 * i6) * i2 * mk
    return [(ff * dx).sum(1), (ff * dy).sum(1), (ff * dz).sum(1)]


def lje(xs, ys, zs):
    """Total potential energy, 4 (r^-12 - r^-6), halved for double counting.

    This is the cross-system check: it is sensitive to every pair distance, so
    two systems agreeing on it are running the same physics.
    """
    dx = xs[:, None] - xs[None, :]
    dy = ys[:, None] - ys[None, :]
    dz = zs[:, None] - zs[None, :]
    r2 = dx * dx + dy * dy + dz * dz + mdbig
    mk = (MDRC2 - r2 >= 0).astype(float)
    i2 = 1.0 / r2
    i6 = i2 * i2 * i2
    return float((4.0 * (i6 * i6 - i6) * mk).sum() / 2.0)


def mdstep(st):
    """One velocity-Verlet step: half-kick, drift, force, half-kick."""
    p, v, f = st
    vh = [v[i] + 0.5 * MDDT * f[i] for i in range(3)]
    pn = [p[i] + MDDT * vh[i] for i in range(3)]
    fn = ljf(pn[0], pn[1], pn[2])
    return [pn, [vh[i] + 0.5 * MDDT * fn[i] for i in range(3)], fn]


def mdrun(m):
    st = nest(mdstep,
              [[mdx, mdy, mdz], [mdzero, mdzero, mdzero], ljf(mdx, mdy, mdz)],
              m)
    return lje(st[0][0], st[0][1], st[0][2])


# Cell-list binning: the irregular half of a real MD code.  100000 atoms are
# assigned to one of 8^3 cells, sorted by cell, and counted -- which is what
# builds the neighbour-list offsets.
NMC = 100000
_j = np.arange(1, NMC + 1, dtype=float)
mcx = 8.0 * np.mod(1.7 * _j, 1.0)
mcy = 8.0 * np.mod(2.3 * _j, 1.0)
mcz = 8.0 * np.mod(3.1 * _j, 1.0)


def mdcells():
    ci = (np.floor(mcx) + 8 * np.floor(mcy) + 64 * np.floor(mcz) + 1).astype(np.int64)
    _, counts = np.unique(np.sort(ci), return_counts=True)
    return int(np.cumsum(counts).sum())


def main():
    print("Experiment 13 -- molecular dynamics")
    print("")

    bench("Lennard-Jones force, 2048 atoms, cut-off",
          lambda: ljf(mdx, mdy, mdz))
    check("Lennard-Jones force", lje(mdx, mdy, mdz))

    bench("velocity-Verlet, 2048 atoms, 10 steps", lambda: mdrun(10))
    check("velocity-Verlet", mdrun(10))

    bench("cell-list binning, 100000 atoms", mdcells)
    check("cell-list binning", mdcells())

    # The .m file's mixed-dtype tuple probe has no NumPy counterpart: a Python
    # list of arrays is a list of references and costs nothing to build,
    # whatever the dtypes.  It is measured here anyway, as the control that
    # shows the Mathilda finding is a representation artefact and not
    # intrinsic to the operation.
    print("")
    print("-- the mixed-dtype tuple return (600 x 600), NumPy control --")
    tta = np.random.rand(600, 600)
    ttb = np.random.rand(600, 600)
    ttm = (tta - 0.5 >= 0).astype(np.int64)
    bench("[a, b]       -- uniform float64", lambda: [tta, ttb])
    bench("[a, b, mask] -- one int64 array", lambda: [tta, ttb, ttm])
    ttu = [tta, ttb]
    bench("caller's next op after the tuple", lambda: ttu[0] + ttu[1])


if __name__ == "__main__":
    main()
