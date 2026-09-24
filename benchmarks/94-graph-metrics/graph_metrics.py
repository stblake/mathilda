#!/usr/bin/env python3
"""Experiment 94 -- Graph distances, centralities and clustering (networkx column).

Same cases as ``graph_metrics.m``, same deterministic graphs, same order.
networkx is pure Python (its PageRank, eigenvector and Katz use scipy/numpy
where networkx does), so -- as in experiments 29 and 92 -- a Mathilda win here
is a weaker claim than a win against a compiled library.

Idioms are the ones a networkx user would reach for, with the options that make
the value match Wolfram's definition: ``betweenness_centrality(normalized=
False)``, ``closeness_centrality(wf_improved=False)``, ``pagerank`` /
``katz_centrality`` / ``eigenvector_centrality_numpy`` (rescaled to total 1),
``all_pairs_shortest_path_length`` for the distance matrix, ``triangles``,
``transitivity`` and ``clustering``. Exact rational checks are computed outside
the timed region from integer triangle counts.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
from fractions import Fraction
import networkx as nx

from harness import check, require

require(["networkx", "networkx:betweenness_centrality", "networkx:closeness_centrality",
         "networkx:pagerank", "networkx:eigenvector_centrality_numpy", "networkx:diameter",
         "networkx:triangles", "networkx:transitivity", "networkx:clustering",
         "networkx:katz_centrality"])


def sw_edges(n):
    e = [(i, i % n + 1) for i in range(1, n + 1)]
    e += [(i, (i + 1 + (7919 * i * i + 104729 * i) % (n // 2 - 2)) % n + 1) for i in range(1, n + 1)]
    return e


nw = 100000
web_edges = [(i, (i + t * 19999 + (7919 * i * i + 104729 * i * t + 31 * t) % 19998) % nw + 1)
             for i in range(1, nw + 1) for t in range(5)]
band_edges = [(i, i + d) for i in range(1, 20001) for d in range(1, 6) if i + d <= 20000]
e10k = sw_edges(10000)
e3k = sw_edges(3000)


def mk(n, edges, directed=False):
    g = nx.DiGraph() if directed else nx.Graph()
    g.add_nodes_from(range(1, n + 1))
    g.add_edges_from(edges)
    return g


new_sw10k = lambda: mk(10000, e10k)
new_sw3k = lambda: mk(3000, e3k)
new_web = lambda: mk(nw, web_edges, True)
new_band = lambda: mk(20000, band_edges)
sw10k, sw3k, web, band = new_sw10k(), new_sw3k(), new_web(), new_band()


def bc(g): return nx.betweenness_centrality(g, normalized=False)
def cc(g): return nx.closeness_centrality(g, wf_improved=False)
def gdm(g): return dict(nx.all_pairs_shortest_path_length(g))
def pr(g): return nx.pagerank(g, alpha=0.85, tol=1e-15, max_iter=10000)
def eig(g):
    v = nx.eigenvector_centrality_numpy(g)
    s = sum(v.values())
    return {k: x / s for k, x in v.items()}
def katz(g): return nx.katz_centrality(g, alpha=0.1, beta=1.0, normalized=False, tol=1e-14, max_iter=10000)
def tri(g): return sum(nx.triangles(g).values()) // 3


def gcc_exact(g):
    t = sum(nx.triangles(g).values())
    trip = sum(d * (d - 1) // 2 for _, d in g.degree())
    f = Fraction(t, trip)
    return "%d/%d" % (f.numerator, f.denominator) if f.denominator != 1 else str(f.numerator)


def lcc_exact(g):
    tr = nx.triangles(g)
    tot = Fraction(0)
    for v, d in g.degree():
        if d >= 2 and tr[v]:
            tot += Fraction(tr[v], d * (d - 1) // 2)
    return "%d/%d" % (tot.numerator, tot.denominator) if tot.denominator != 1 else str(tot.numerator)


def timed(f, g, reps):
    """Minimum wall time over `reps` runs (after one untimed warm-up when
    reps > 1, as harness.bench does) and the result of the last run."""
    if reps > 1:
        f(g)
    best, out = None, None
    for _ in range(reps):
        t0 = time.perf_counter()
        out = f(g)
        dt = time.perf_counter() - t0
        best = dt if best is None or dt < best else best
    return best, out


# (label, fresh-graph maker, graph, timed function, check-from-result, reps).
# The O(nm) cases take minutes per run in networkx, so they are timed once and
# their check is read off the timed result rather than recomputed.
CASES = [
    ("BetweennessCentrality, 10^4-vertex small world", new_sw10k, sw10k, bc,
     lambda r, g: str(round(sum(r.values()))), 1),
    ("ClosenessCentrality, 10^4-vertex small world", new_sw10k, sw10k, cc,
     lambda r, g: str(round(10**6 * sum(r.values()))), 1),
    ("GraphDiameter, 10^4-vertex small world", new_sw10k, sw10k, nx.diameter,
     lambda r, g: str(r), 1),
    ("GraphDistanceMatrix, 3000-vertex small world", new_sw3k, sw3k, gdm,
     lambda r, g: str(sum(sum(x.values()) for x in r.values())), 1),
    ("PageRankCentrality, 10^5-vertex web graph", new_web, web, pr,
     lambda r, g: str(round(10**8 * max(r.values()))), 3),
    ("EigenvectorCentrality, 10^4-vertex small world", new_sw10k, sw10k, eig,
     lambda r, g: str(round(10**8 * max(r.values()))), 3),
    ("KatzCentrality, 10^5-vertex web graph", new_web, web, katz,
     lambda r, g: str(round(10**9 * max(r.values()))), 3),
    ("GraphTriangleCount, 10^5-edge band graph", new_band, band, tri,
     lambda r, g: str(r), 3),
    ("GlobalClusteringCoefficient, 10^5-edge band graph", new_band, band, nx.transitivity,
     lambda r, g: gcc_exact(g), 3),
    ("LocalClusteringCoefficient, 10^5-edge band graph", new_band, band, nx.clustering,
     lambda r, g: lcc_exact(g), 3),
]

for label, mkg, g, f, chk, reps in CASES:
    t, r = timed(f, g, reps)
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * t), flush=True)
    check(label, chk(r, g))

# ---- COLD cases: a fresh graph object per timed call (construction untimed).
# networkx caches nothing on a graph, so these re-measure the warm cost; they
# are kept so the three columns line up case for case.
for label, mkg, g, f, chk, reps in CASES:
    best, r, gg = None, None, None
    for _ in range(reps):
        gg = mkg()
        t0 = time.perf_counter()
        r = f(gg)
        dt = time.perf_counter() - t0
        best = dt if best is None or dt < best else best
    print("BENCH\t%s (cold)\t%.3f" % (label, 1000.0 * best), flush=True)
    check(label + " (cold)", chk(r, gg))
