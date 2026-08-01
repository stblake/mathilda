#!/usr/bin/env python3
"""Experiment 19 -- Ray tracing: branch-free arrays, and a 65-element table.

Runs the same algorithm as ``ray_tracing.m``, in the same order, with the same
sizes.  See ``README.md`` for the measurements and the analysis.

    python3 ray_tracing.py

WHAT IT MEASURES.  A ray tracer written for an array language is BRANCH-FREE
by construction: there is no ``if`` anywhere, because every conditional
becomes a mask and every selection becomes arithmetic on that mask.  That
style is how all high-throughput array code expresses control flow.

It also ends with the operation experiment 12 is about: once the winning
sphere per ray is known, its centre must be GATHERED by index -- from a
65-entry table, with 262144 indices.  That is the shape no packing threshold
can help with, because the table is small *by design*.

DETERMINISM.  The scene, the camera and the light are all closed forms, so the
mean pixel intensity is an exact cross-system check.
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


# ---- camera and scene ----------------------------------------------------

RTW, RTH = 512, 512
RTNP = RTW * RTH           # 262144 rays, one per pixel
RTNS = 64                  # spheres
RTR = 0.35
RTR2 = RTR ** 2

# Every ray is flattened into a plain vector, so the whole image is one array
# operation per sphere rather than a loop over pixels.
_ri = np.arange(1, RTH + 1)[:, None]
_rj = np.arange(1, RTW + 1)[None, :]
rtjs = np.broadcast_to((_rj - (RTW + 1) / 2.0) / RTW, (RTH, RTW)).ravel().copy()
rtis = np.broadcast_to((_ri - (RTH + 1) / 2.0) / RTH, (RTH, RTW)).ravel().copy()
rtnm = 1.0 / np.sqrt(rtjs ** 2 + rtis ** 2 + 1.0)
rtdx = rtjs * rtnm
rtdy = rtis * rtnm
rtdz = rtnm

_s = np.arange(1, RTNS + 1)
rtcx = (-3.5 + (_s - 1) % 8).astype(float)
rtcy = (-3.5 + (_s - 1) // 8).astype(float)
rtcz = 6.0 + 0.4 * np.sin(1.7 * _s)

# A dummy entry at position 0 so a MISS (sphere id 0) can index the table
# without a branch: cid is then always in range, and the miss is masked out at
# the end instead.
rtcxp = np.concatenate(([0.0], rtcx))
rtcyp = np.concatenate(([0.0], rtcy))
rtczp = np.concatenate(([0.0], rtcz))

RTLX, RTLY, RTLZ = 0.577, -0.577, -0.577       # light direction


def rtrender():
    tmin = np.full(RTNP, 1e9)                  # nearest hit so far, per ray
    cid = np.zeros(RTNP, dtype=np.int64)       # which sphere; 0 means a miss
    s = 1
    while s <= RTNS:
        cx, cy, cz = rtcx[s - 1], rtcy[s - 1], rtcz[s - 1]

        # Ray-sphere intersection, from the eye at the origin.  tca is the
        # projection of the centre onto the ray; disc < 0 means a miss.
        tca = rtdx * cx + rtdy * cy + rtdz * cz
        d2 = (cx * cx + cy * cy + cz * cz) - tca * tca
        disc = RTR2 - d2

        # Clip before the sqrt so a miss produces 0 rather than a NaN.
        tt = tca - np.sqrt(np.clip(disc, 0.0, 1e9))

        ok = (disc >= 0).astype(np.int64) * (tt - 0.001 >= 0).astype(np.int64)

        # A miss is pushed to +infinity so the running minimum ignores it.
        cand = ok * tt + (1 - ok) * 1e9

        # STRICTLY closer: the epsilon matters, because at s == 1 a missing
        # ray has cand == tmin == 1e9 exactly, and a non-strict test would
        # make every missed ray claim to have hit sphere 1.
        bt = (tmin - cand - 1e-6 >= 0).astype(np.int64)
        cid = bt * s + (1 - bt) * cid          # select, no branch
        tmin = np.minimum(tmin, cand)
        s += 1

    hx, hy, hz = rtdx * tmin, rtdy * tmin, rtdz * tmin      # hit points

    # THE GATHER: 262144 indices into a 65-entry table.
    gx, gy, gz = rtcxp[cid], rtcyp[cid], rtczp[cid]

    nx, ny, nz = (hx - gx) / RTR, (hy - gy) / RTR, (hz - gz) / RTR   # normals
    lam = np.clip(nx * RTLX + ny * RTLY + nz * RTLZ, 0.0, 1.0)       # Lambert
    vis = (1e8 - tmin >= 0).astype(float)                            # hit mask
    return float((lam * vis).sum() / RTNP)


def main():
    print("Experiment 19 -- ray tracing")
    print("")

    bench("ray trace, 512^2 rays x 64 spheres, diffuse", rtrender)
    check("ray trace (mean pixel intensity)", rtrender())

    # The gather, alone.  The source table has 65 entries -- small BY DESIGN,
    # so no packing threshold can help -- and the index array has 262144.
    print("")
    print("-- the sphere-parameter lookup, alone --")
    rtcid = np.random.randint(0, 65, RTNP)
    bench("table[cid], 262144 indices from 65", lambda: rtcxp[rtcid])

    # One sphere's worth of the inner loop, to show where the rest goes:
    # fifteen whole-array passes over 262144 float64 is ~32 MB of traffic per
    # sphere, and there are 64 spheres.
    print("")
    print("-- one sphere of the inner loop --")
    rttca = rtdx * 1.0 + rtdy * 2.0 + rtdz * 6.0
    bench("tca = dx cx + dy cy + dz cz",
          lambda: rtdx * 1.0 + rtdy * 2.0 + rtdz * 6.0)
    bench("sqrt(clip(disc, ...))",
          lambda: np.sqrt(np.clip(RTR2 - rttca, 0.0, 1e9)))
    bench("(disc >= 0) * (tt >= 0)",
          lambda: (rttca >= 0).astype(np.int64) * (rttca - 0.001 >= 0).astype(np.int64))
    bench("minimum(a, b)", lambda: np.minimum(rttca, rtdz))


if __name__ == "__main__":
    main()
