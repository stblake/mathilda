#!/usr/bin/env python3
"""Experiment 96 -- the native hypergraph subsystem (Python column).

Same cases, same deterministic hypergraphs and same labels as
``hypergraphs.m``.  The library is xgi 0.10 (``pip install xgi``; run with the
interpreter that has it, e.g. a venv), a pure-Python hypergraph package over
networkx/numpy/scipy.  Where xgi has no such function the idiom is the obvious
one on top of it (line graph + networkx BFS for s-walk distances, a DiGraph
for the FR HypergraphToGraph).  The minimum hitting set is solved as a 0-1
integer program with scipy's HiGHS (``scipy.optimize.milp``), the strongest
off-the-shelf solver a Python user has.  Neither xgi nor networkx enumerates
minimal transversals, so that case is SKIPped (absent), not approximated.

xgi builds hyperedges as SETS, matching the set semantics every set-based
Mathilda head uses.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import random
import time

from harness import bench, bench_if, check, require

require(["xgi", "networkx", "numpy", "scipy:optimize"])

try:
    import xgi
    import networkx as nx
    import numpy as np
except ImportError:          # the SKIP lines below still name every case
    xgi = None


def mk_e(n, m, p, q):
    out = []
    for i in range(1, m + 1):
        a = (7919 * i) % n
        b = (a + 1 + i % p) % n
        c = (b + 1 + (3 * i) % q) % n
        out.append([a + 1, b + 1, c + 1])
    return out


def lcg(n, m):
    x, vals = 42, []
    for _ in range(3 * m):
        x = (1103515245 * x + 12345) % 2147483648
        vals.append((x // 65536) % n + 1)
    return [vals[k:k + 3] for k in range(0, len(vals), 3)]


nB, eB = 20000, mk_e(20000, 100000, 101, 103)
nL, eL = 100000, mk_e(100000, 100000, 101, 103)
nI, eI = 1000, mk_e(1000, 10000, 101, 103)
eM = lcg(80, 200)


def build(n, e):
    h = xgi.Hypergraph()
    h.add_nodes_from(range(1, n + 1))
    h.add_edges_from(e)
    return h


def wl(b):
    return "True" if b else "False"


def wl_list(xs):
    return "{" + ", ".join(str(x) for x in xs) + "}"


def finite(d):
    return wl_list([max(d.values()), len(d)])


def to_graph_fr(e):
    g = nx.DiGraph()
    for s in e:
        for a in range(len(s)):
            for b in range(a + 1, len(s)):
                if s[a] != s[b]:
                    g.add_edge(s[a], s[b])
    return g


def edge_distances(h, s=1):
    lg = xgi.to_line_graph(h, s=s)
    first = min(h.edges)
    return nx.single_source_shortest_path_length(lg, first)


def vertex_distances(h):
    # xgi.single_source_shortest_path_length is an O(n^2) Dijkstra (linear
    # min-scan) -- hours at n = 10^5 -- so the idiom is BFS on the clique
    # expansion, exactly what the Mathematica column does.
    return nx.single_source_shortest_path_length(xgi.to_graph(h), 1)


def min_transversal(e):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    verts = sorted({v for s in e for v in s})
    pos = {v: i for i, v in enumerate(verts)}
    a = lil_matrix((len(e), len(verts)))
    for r, s in enumerate(e):
        for v in set(s):
            a[r, pos[v]] = 1
    res = milp(np.ones(len(verts)), constraints=LinearConstraint(a.tocsr(), lb=1),
               integrality=np.ones(len(verts)), bounds=Bounds(0, 1))
    return [v for v, x in zip(verts, res.x) if x > 0.5]


if xgi is not None:
    hB, hL, hI = build(nB, eB), build(nL, eL), build(nI, eI)

X = "xgi"

bench_if("Hypergraph construction, 10^5 hyperedges", X, lambda: build(nB, eB))
check_v = (lambda: build(nB, eB).num_nodes) if xgi else None
if xgi: check("Hypergraph construction, 10^5 hyperedges", check_v())

bench_if("VertexDegree, 10^5 hyperedges", X, lambda: hB.nodes.degree.asnumpy())
if xgi: check("VertexDegree, 10^5 hyperedges", int(hB.nodes.degree.asnumpy().sum()))

bench_if("HyperedgeSizes, 10^5 hyperedges", X, lambda: hB.edges.size.asnumpy())
if xgi: check("HyperedgeSizes, 10^5 hyperedges", int(hB.edges.size.asnumpy().sum()))

bench_if("HypergraphDual, 10^5 hyperedges", X, lambda: hB.dual())
if xgi: check("HypergraphDual, 10^5 hyperedges", int(hB.dual().edges.size.asnumpy().sum()))

bench_if("HypergraphCliqueExpansion, 10^5 hyperedges", X, lambda: xgi.to_graph(hB))
if xgi: check("HypergraphCliqueExpansion, 10^5 hyperedges", xgi.to_graph(hB).number_of_edges())

bench_if("HypergraphStarExpansion, 10^5 hyperedges", X, lambda: xgi.to_bipartite_graph(hB))
if xgi: check("HypergraphStarExpansion, 10^5 hyperedges", xgi.to_bipartite_graph(hB).number_of_edges())

bench_if("HypergraphToGraph (FR), 10^5 hyperedges", X, lambda: to_graph_fr(eB))
if xgi: check("HypergraphToGraph (FR), 10^5 hyperedges", to_graph_fr(eB).number_of_edges())

bench_if("HypergraphLineGraph, 10^5 hyperedges", X, lambda: xgi.to_line_graph(hL))
if xgi: check("HypergraphLineGraph, 10^5 hyperedges", xgi.to_line_graph(hL).number_of_edges())

bench_if("HypergraphConnectedComponents, 10^5 hyperedges", X, lambda: list(xgi.connected_components(hL)))
if xgi: check("HypergraphConnectedComponents, 10^5 hyperedges",
              wl_list(sorted(len(c) for c in xgi.connected_components(hL))))

bench_if("ConnectedHypergraphQ (FR), 10^5 hyperedges", X, lambda: xgi.is_connected(build(nB, eB)))
if xgi: check("ConnectedHypergraphQ (FR), 10^5 hyperedges", wl(xgi.is_connected(hB)))

bench_if("HyperedgeConnectedComponents s=2, 10^5 hyperedges", X,
         lambda: nx.number_connected_components(xgi.to_line_graph(hB, s=2)))
if xgi: check("HyperedgeConnectedComponents s=2, 10^5 hyperedges",
              nx.number_connected_components(xgi.to_line_graph(hB, s=2)))

bench_if("HyperedgeDistance from e1, 10^5 hyperedges", X, lambda: edge_distances(hL))
if xgi: check("HyperedgeDistance from e1, 10^5 hyperedges", finite(edge_distances(hL)))

bench_if("HypergraphDistance from v1, 10^5 hyperedges", X,
         lambda: vertex_distances(hL))
if xgi: check("HypergraphDistance from v1, 10^5 hyperedges",
              finite(vertex_distances(hL)))

bench_if("IncidenceMatrix, 1000 x 10^4", X, lambda: xgi.incidence_matrix(hI, sparse=False))
if xgi: check("IncidenceMatrix, 1000 x 10^4", int(xgi.incidence_matrix(hI, sparse=False).sum()))


def fr_random():
    r = random.Random(1)
    return build(nB, [[r.randint(1, nB) for _ in range(3)] for _ in range(100000)])


bench_if("RandomHypergraph (FR form), 10^5 triples", X, fr_random)
if xgi: check("RandomHypergraph (FR form), 10^5 triples", fr_random().num_edges)

print("SKIP\tTransversalHypergraph (FR), 16 vertices\txgi:transversal", flush=True)

bench_if("FindMinimumTransversal, 80 vertices x 200", "scipy:optimize", lambda: min_transversal(eM))
check("FindMinimumTransversal, 80 vertices x 200", len(min_transversal(eM)))


# ---- COLD cases (see hypergraphs.m): fresh object per timed call ----------
def cold_bench(label, mk, f, reps=3):
    if xgi is None:
        print("SKIP\t%s\txgi" % label, flush=True)
        return
    ts = []
    for _ in range(reps):
        hh = mk()
        t0 = time.perf_counter()
        f(hh)
        ts.append(time.perf_counter() - t0)
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * min(ts)), flush=True)


newB = lambda: build(nB, eB)
newL = lambda: build(nL, eL)
newI = lambda: build(nI, eI)

cold_bench("VertexDegree, 10^5 hyperedges (cold)", newB, lambda h: h.nodes.degree.asnumpy())
if xgi: check("VertexDegree, 10^5 hyperedges (cold)", int(newB().nodes.degree.asnumpy().sum()))
cold_bench("HyperedgeSizes, 10^5 hyperedges (cold)", newB, lambda h: h.edges.size.asnumpy())
if xgi: check("HyperedgeSizes, 10^5 hyperedges (cold)", int(newB().edges.size.asnumpy().sum()))
cold_bench("HypergraphDual, 10^5 hyperedges (cold)", newB, lambda h: h.dual())
if xgi: check("HypergraphDual, 10^5 hyperedges (cold)", int(newB().dual().edges.size.asnumpy().sum()))
cold_bench("HypergraphCliqueExpansion, 10^5 hyperedges (cold)", newB, xgi.to_graph if xgi else None)
if xgi: check("HypergraphCliqueExpansion, 10^5 hyperedges (cold)", xgi.to_graph(newB()).number_of_edges())
cold_bench("HypergraphStarExpansion, 10^5 hyperedges (cold)", newB, xgi.to_bipartite_graph if xgi else None)
if xgi: check("HypergraphStarExpansion, 10^5 hyperedges (cold)",
              xgi.to_bipartite_graph(newB()).number_of_edges())
cold_bench("HypergraphLineGraph, 10^5 hyperedges (cold)", newL, lambda h: xgi.to_line_graph(h))
if xgi: check("HypergraphLineGraph, 10^5 hyperedges (cold)", xgi.to_line_graph(newL()).number_of_edges())
cold_bench("HypergraphConnectedComponents, 10^5 hyperedges (cold)", newL,
           lambda h: list(xgi.connected_components(h)))
if xgi: check("HypergraphConnectedComponents, 10^5 hyperedges (cold)",
              wl_list(sorted(len(c) for c in xgi.connected_components(newL()))))
cold_bench("HyperedgeConnectedComponents s=2, 10^5 hyperedges (cold)", newB,
           lambda h: nx.number_connected_components(xgi.to_line_graph(h, s=2)))
if xgi: check("HyperedgeConnectedComponents s=2, 10^5 hyperedges (cold)",
              nx.number_connected_components(xgi.to_line_graph(newB(), s=2)))
cold_bench("HyperedgeDistance from e1, 10^5 hyperedges (cold)", newL, edge_distances)
if xgi: check("HyperedgeDistance from e1, 10^5 hyperedges (cold)", finite(edge_distances(newL())))
cold_bench("HypergraphDistance from v1, 10^5 hyperedges (cold)", newL,
           vertex_distances)
if xgi: check("HypergraphDistance from v1, 10^5 hyperedges (cold)",
              finite(vertex_distances(newL())))
cold_bench("IncidenceMatrix, 1000 x 10^4 (cold)", newI, lambda h: xgi.incidence_matrix(h, sparse=False))
if xgi: check("IncidenceMatrix, 1000 x 10^4 (cold)", int(xgi.incidence_matrix(newI(), sparse=False).sum()))
