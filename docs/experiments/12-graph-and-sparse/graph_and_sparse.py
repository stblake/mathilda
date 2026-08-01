#!/usr/bin/env python3
"""Experiment 12 -- Graph analytics and sparse linear algebra (NumPy column).

Runs the same algorithm as ``graph_and_sparse.m``, in the same order, with the
same sizes.  See ``README.md`` for the measurements and the analysis.

    python3 graph_and_sparse.py

WHAT IT MEASURES.  Every sparse and graph kernel is built out of one
operation: gather -- ``y = x[idx]``, reading a value array through an index
array.  A sparse matrix-vector product is a gather, a multiply and a segmented
sum; a PageRank iteration is twenty of those; a breadth-first search is a
gather of adjacency rows and a set difference.

REPRESENTATION.  The graph is a dense ``gn x gdeg`` matrix of neighbour ids --
CSR with the row pointers implied, which is what a uniform-degree graph gives
you.  NumPy is given exactly the same representation as the two CAS rather
than ``scipy.sparse``, because the point of the row is to compare the *gather*,
not to compare three sparse-matrix libraries.

INDEXING.  The Wolfram Language is 1-based and NumPy is 0-based, so every
index array built to match the ``.m`` file is used as ``x[idx - 1]``.  Keeping
the arrays 1-based (rather than generating 0-based ones here) means the two
files hold the same numbers, which matters when a check value disagrees and
has to be traced.

DETERMINISM.  The timed graph is random; the three systems cannot be made to
draw the same one, so every check runs a separate deterministic 4096-node
graph.  A timing is meaningless until the three systems agree on the check.
"""

import time

import numpy as np

# ---- shared reporting helpers (identical in every experiment file) --------


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs.

    The minimum, not the mean: we are measuring the cost of the work, and
    every source of noise on a loaded machine can only add.
    """
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


# ---- the timed graph -----------------------------------------------------

GN = 100000                            # nodes
GDEG = 16                              # in-degree, uniform

# gadj[i] is the list of nodes that link TO node i.  Uniform in-degree is what
# makes the row pointers implicit: row i occupies flat positions
# i*GDEG .. (i+1)*GDEG.
gadj = np.random.randint(1, GN + 1, (GN, GDEG))
gflat = gadj.ravel()                   # the CSR column-index array
gv = np.random.rand(GN)                # a value per node
gw = np.random.rand(GN * GDEG)         # a weight per edge

# ---- the deterministic check graph ---------------------------------------

GCN = 4096
GCD = 8
gca = (np.arange(1, GCN + 1)[:, None] * np.arange(1, GCD + 1)[None, :]
       + np.arange(1, GCD + 1)[None, :]) % GCN + 1
gcf = gca.ravel()
gcw = ((np.arange(1, GCN * GCD + 1) % 7) + 1) / 7.0


# ---- kernels -------------------------------------------------------------


def pagerank(f, nn, dd, m):
    """PageRank by power iteration.

        p[i] <- (1-d)/n + d * sum over the in-neighbours j of i of p[j]/deg

    Three whole-array operations, which is the point: the gather ``p[f-1]``
    produces one value per EDGE, ``reshape`` groups those into one row per
    node, and ``sum(1)`` sums each row.  No scalar loop and no sparse-matrix
    type -- just the gather and a segmented reduction.
    """
    p = np.full(nn, 1.0 / nn)
    for _ in range(m):
        p = 0.15 / nn + 0.85 * p[f - 1].reshape(nn, dd).sum(1) / dd
    return p


def bfs(a, start, levels):
    """Breadth-first search: level-synchronous frontier expansion.

    ``a[fr - 1]`` gathers whole adjacency ROWS -- a rank-2 gather -- and the
    set operations do the work of a visited-bitmap.  This is the Graph500
    kernel written the way an array language expresses it.
    """
    visited = np.array([start])
    frontier = np.array([start])
    k = 0
    while k < levels and frontier.size > 0:
        nxt = np.setdiff1d(np.unique(a[frontier - 1].ravel()), visited)
        visited = np.union1d(visited, nxt)
        frontier = nxt
        k += 1
    return int(visited.size)


def main():
    print("Experiment 12 -- graph analytics and sparse linear algebra")
    print("")

    # 1. The primitive, alone: 1.6e6 indices into a 100000-element vector.
    bench("gather, 1.6e6 indices into 100000", lambda: gv[gflat - 1])
    check("gather", float(gcw[(7 * np.arange(1, 32769)) % 4096].sum()))

    # 2. One sparse matrix-vector product: gather, multiply, segmented sum.
    bench("SpMV, 100000 x 16 CSR",
          lambda: (gw * gv[gflat - 1]).reshape(GN, GDEG).sum(1))
    check("SpMV",
          float((gcw * gcw[gcf - 1]).reshape(GCN, GCD).sum(1).sum()))

    # 3. Twenty of the above, with the damping term.
    bench("PageRank, 1.6e6 edges, 20 iterations",
          lambda: pagerank(gflat, GN, GDEG, 20))
    check("PageRank", float(pagerank(gcf, GCN, GCD, 20).sum()))

    # 4. Five levels of frontier expansion from one source.
    bench("breadth-first search, 5 levels", lambda: bfs(gadj, 1, 5))
    check("breadth-first search", bfs(gca, 1, 5))


if __name__ == "__main__":
    main()
