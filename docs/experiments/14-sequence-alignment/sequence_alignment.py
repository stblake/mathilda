#!/usr/bin/env python3
"""Experiment 14 -- Bioinformatics: sequence alignment and k-mer counting.

Runs the same algorithm as ``sequence_alignment.m``, in the same order, with
the same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 sequence_alignment.py

WHAT IT MEASURES.  Everything here is INTEGER end to end -- bases, scores, gap
penalties, k-mer keys, and an alignment matrix of integer dynamic programming.
That is what makes it a useful probe of an array system: if the integer arms
of the elementwise, scan and structural paths are missing, this finds them and
a floating-point workload never will.

THE ALIGNMENT, AND WHY IT IS WRITTEN THIS WAY.  The Needleman-Wunsch
within-row gap recurrence

    h[j] = max(cand[j], h[j-1] - g)

is a max-plus prefix scan and looks sequential.  Substituting
``H[j] = h[j] + g*j`` turns it into a plain running maximum, so all three
systems run *that* -- ``np.maximum.accumulate`` here, ``FoldList[Max, ...]``
in the two CAS.  The columns then compare execution rather than cleverness.

DETERMINISM.  Both sequences and the k-mer input come from closed forms, so
every check is an integer that all three systems must agree on exactly.

WHERE NUMPY HAS AN UNFAIR ADVANTAGE, STATED PLAINLY.
``sliding_window_view`` is a stride trick: it returns a *view* and copies
nothing, where both CAS materialise 6 million elements.  That is a real
capability difference, not a tuning difference, and it is why the k-mer row
should be read as "NumPy has strided views and Mathilda does not" rather than
as "Mathilda's Partition is slow".
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


# ---- 1. Needleman-Wunsch global alignment --------------------------------

BM = 2000                  # sequence length; the DP matrix is BM x BM
BGAP = 1                   # linear gap penalty

_b = np.arange(1, BM + 1, dtype=np.int64)
bs1 = (_b * _b + 3 * _b) % 4                    # bases 0..3
bs2 = (7 * _b + _b // 3) % 4

brng1 = np.arange(1, BM + 1, dtype=np.int64)    # 1..BM, the g*j shift
brng0 = np.arange(0, BM + 1, dtype=np.int64)    # 0..BM, the g*j unshift


def nwrow(prev, ch):
    """One row of the DP matrix, from the row above it.

    ``sco`` is the substitution score, +2 on a match and -1 on a mismatch,
    computed without a comparison predicate: for integer bases
    ``abs(b - ch)`` is 0 exactly when they are equal, so
    ``(abs(b - ch) - 1 >= 0)`` is the "not equal" indicator.  Written with a
    Python-level conditional it would be one call per column.
    """
    sco = 2 - 3 * (np.abs(bs2 - ch) - 1 >= 0).astype(np.int64)
    dg = prev[:-1] + sco                        # diagonal: match/mismatch
    up = prev[1:] - BGAP                        # from above: a gap
    cand = np.maximum(dg, up) + BGAP * brng1
    hh = np.maximum.accumulate(
        np.concatenate(([prev[0] - BGAP], cand)))   # the running maximum
    return hh - BGAP * brng0                    # undo the shift


def nwalign():
    row = -BGAP * brng0
    for ch in bs1:
        row = nwrow(row, ch)
    return int(row[-1])


# ---- 2. k-mer encoding and distinct count --------------------------------

BKN = 500000               # bases
BK = 12                    # k-mer length; 4^12 = 16.7e6 possible keys

_k = np.arange(1, BKN + 1, dtype=np.int64)
bcode = (_k * _k + 3 * _k + _k // 7) % 4
bpow = (4 ** (BK - np.arange(1, BK + 1))).astype(np.int64)   # base-4 places


def kmers():
    """Sliding window, base-4 encode, distinct count.

    ``sliding_window_view`` returns a view over ``bcode`` -- no copy at all --
    which is the whole of NumPy's advantage on this row.
    """
    w = np.lib.stride_tricks.sliding_window_view(bcode, BK)
    return int(np.unique(w @ bpow).size)


# ---- 3. rolling GC content -----------------------------------------------

BGN = 10 ** 7
BGW = 1000                 # window
bgseq = np.random.randint(0, 2, BGN)            # timed input: random
_g = np.arange(1, 100001, dtype=np.int64)
bgcs = (_g ** 3 + _g + _g // 3) % 2             # check input: deterministic


def gcwin(s, w):
    """A rolling sum is two shifted reads of the prefix sum -- the standard
    trick, and the reason this row is a scan benchmark and not a convolution
    one."""
    cs = np.cumsum(s)
    return int((cs[w:] - cs[:-w]).sum())


def main():
    print("Experiment 14 -- sequence alignment and k-mer counting")
    print("")

    bench("Needleman-Wunsch alignment, 2000 x 2000", nwalign)
    check("Needleman-Wunsch alignment", nwalign())

    bench("k-mer encode + distinct count, 5e5, k=12", kmers)
    check("k-mer encode + distinct count", kmers())

    bench("rolling GC content, 10^7 bases, w=1000", lambda: gcwin(bgseq, BGW))
    check("rolling GC content", gcwin(bgcs, 100))

    # The primitives underneath, for the same split the .m file reports.
    print("")
    print("-- the inner row's primitives, per 2000-element call --")
    nwprev = -np.arange(0, BM + 1, dtype=np.int64)

    def rep(f):
        return lambda: [f() for _ in range(200)]

    bench("abs(seq - ch)", rep(lambda: np.abs(bs2 - 2)))
    bench("(x >= 0) mask", rep(lambda: (bs2 - 2 >= 0).astype(np.int64)))
    bench("prev[:-1]", rep(lambda: nwprev[:-1].copy()))
    bench("maximum(prev[:-1], prev[1:])",
          rep(lambda: np.maximum(nwprev[:-1], nwprev[1:])))
    bench("maximum.accumulate(v)", rep(lambda: np.maximum.accumulate(brng1)))
    bench("the whole row", rep(lambda: nwrow(nwprev, 2)))

    print("")
    print("-- the k-mer pipeline, split --")
    bench("sliding_window_view (a VIEW, no copy)",
          lambda: np.lib.stride_tricks.sliding_window_view(bcode, BK))
    w = np.lib.stride_tricks.sliding_window_view(bcode, BK)
    bench("... @ pow4 (the encode)", lambda: w @ bpow)
    keys = w @ bpow
    bench("unique (the distinct count)", lambda: np.unique(keys))


if __name__ == "__main__":
    main()
