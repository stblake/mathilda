#!/usr/bin/env python3
"""Experiment 93 -- Graph editing, transforms, set operations (networkx column).

Same cases as ``graph_ops_editing.m``, same graphs, same order.  Mathematica and
Mathilda graphs are immutable values, so every mutating networkx idiom is timed
together with the ``copy()`` that makes it produce a new graph -- the three
columns all time "produce the edited graph".  networkx is pure Python, so a
Mathilda win here is a weaker claim than one against a compiled library.

Idioms: ``subgraph(...).copy()`` for Subgraph, ``compose`` for GraphUnion
(shared vertices), ``intersection`` / ``difference`` (same vertex sets here),
``disjoint_union``, ``complement``, ``reverse``, ``to_undirected`` /
``to_directed``, ``line_graph``, ``eulerian_circuit``, ``shortest_path`` for
FindPath (any s-t path), ``find_cycle``, ``convert_node_labels_to_integers``.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
import networkx as nx

from harness import bench, check, require

require(["networkx", "networkx:compose", "networkx:intersection", "networkx:difference",
         "networkx:disjoint_union", "networkx:complement", "networkx:line_graph",
         "networkx:eulerian_circuit", "networkx:find_cycle"])

nv = 100000
gE = [(i, i + 1) for i in range(1, nv)] + [(i, i + 2 + (7 * i) % 50) for i in range(1, nv - 51)]
g2E = [(i, i + 3) for i in range(1, nv - 2)]
ne = 50000
euE = [(i, i % ne + 1) for i in range(1, ne + 1)] + [(i, (i + 1) % ne + 1) for i in range(1, ne + 1)]


def mk(edges, n, directed=False):
    h = nx.DiGraph() if directed else nx.Graph()
    h.add_nodes_from(range(1, n + 1))
    h.add_edges_from(edges)
    return h


g = mk(gE, nv)
g2 = mk(g2E, nv)
dg = mk(gE, nv, directed=True)
eu = mk(euE, ne)
c1k = nx.cycle_graph(range(1, 1001))
ns = 2000
hE = [(i, i + 1) for i in range(1, ns)] + [(i, i + 2 + (7 * i) % 50) for i in range(1, ns - 51)]
h2E = [(i, i + 3) for i in range(1, ns - 2)]
h = mk(hE, ns)
h2 = mk(h2E, ns)
delV = list(range(1, nv + 1, 100))
addV = list(range(nv + 1, nv + 10001))
addE = [(i, i + 1000) for i in range(1, 10001)]
delE = [(i, i + 1) for i in range(1, 20001, 2)]
subV = list(range(1, 50001))


def cnt(h):
    return "{%d, %d}" % (h.number_of_nodes(), h.number_of_edges())


def vdel(h):
    h = h.copy(); h.remove_nodes_from(delV); return h


def vadd(h):
    h = h.copy(); h.add_nodes_from(addV); return h


def eadd(h):
    h = h.copy(); h.add_edges_from(addE); return h


def edel(h):
    h = h.copy(); h.remove_edges_from(delE); return h


def tour_ok(c, gr):
    a = sorted(tuple(sorted(e)) for e in c)
    b = sorted(tuple(sorted(e)) for e in gr.edges())
    return "{%d, %s}" % (len(c), "True" if a == b else "False")


bench("VertexDelete, 1000 of 100000 vertices", lambda: vdel(g))
check("VertexDelete, 1000 of 100000 vertices", cnt(vdel(g)))
bench("VertexAdd, 10000 vertices", lambda: vadd(g))
check("VertexAdd, 10000 vertices", cnt(vadd(g)))
bench("EdgeAdd, 10000 edges", lambda: eadd(g))
check("EdgeAdd, 10000 edges", cnt(eadd(g)))
bench("EdgeDelete, 10000 edges", lambda: edel(g))
check("EdgeDelete, 10000 edges", cnt(edel(g)))
bench("Subgraph, 50000 of 100000 vertices", lambda: g.subgraph(subV).copy())
check("Subgraph, 50000 of 100000 vertices", cnt(g.subgraph(subV).copy()))
bench("IndexGraph, 100000 vertices", lambda: nx.convert_node_labels_to_integers(g))
check("IndexGraph, 100000 vertices", cnt(nx.convert_node_labels_to_integers(g)))
bench("GraphUnion, 2 x 100000 vertices", lambda: nx.compose(g, g2))
check("GraphUnion, 2 x 100000 vertices", cnt(nx.compose(g, g2)))
bench("GraphIntersection, 2 x 2000 vertices", lambda: nx.intersection(h, h2))
check("GraphIntersection, 2 x 2000 vertices", cnt(nx.intersection(h, h2)))
bench("GraphDifference, 2 x 2000 vertices", lambda: nx.difference(h, h2))
check("GraphDifference, 2 x 2000 vertices", cnt(nx.difference(h, h2)))
bench("GraphDisjointUnion, 2 x 100000 vertices", lambda: nx.disjoint_union(g, g2))
check("GraphDisjointUnion, 2 x 100000 vertices", cnt(nx.disjoint_union(g, g2)))
bench("GraphComplement, CycleGraph[1000]", lambda: nx.complement(c1k))
check("GraphComplement, CycleGraph[1000]", cnt(nx.complement(c1k)))
bench("ReverseGraph, 199947 arcs", lambda: dg.reverse(copy=True))
check("ReverseGraph, 199947 arcs", cnt(dg.reverse(copy=True)))
bench("UndirectedGraph, 199947 arcs", lambda: dg.to_undirected())
check("UndirectedGraph, 199947 arcs", cnt(dg.to_undirected()))
bench("DirectedGraph, 199947 edges", lambda: g.to_directed())
check("DirectedGraph, 199947 edges", cnt(g.to_directed()))
bench("LineGraph, 99997-edge graph", lambda: nx.line_graph(g2))
check("LineGraph, 99997-edge graph", cnt(nx.line_graph(g2)))
bench("FindEulerianCycle, 100000 edges", lambda: list(nx.eulerian_circuit(eu)))
check("FindEulerianCycle, 100000 edges", tour_ok(list(nx.eulerian_circuit(eu)), eu))
bench("FindPath, 100000 vertices", lambda: nx.shortest_path(g, 1, nv))
check("FindPath, 100000 vertices", "{1, %d}" % nx.shortest_path(g, 1, nv)[-1])
bench("FindCycle, 100000 vertices", lambda: nx.find_cycle(g))
check("FindCycle, 100000 vertices", 1 if nx.find_cycle(g) else 0)


# ---- COLD cases (see graph_ops_editing.m) ---------------------------------
# networkx keeps no per-graph cache, but the fresh-graph protocol is mirrored
# exactly so the three columns time the same thing.
def cold_bench(label, mkf, f, reps=5):
    ts = []
    for _ in range(reps):
        gg = mkf()
        t0 = time.perf_counter()
        f(gg)
        ts.append(time.perf_counter() - t0)
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * min(ts)), flush=True)


new_g = lambda: mk(gE, nv)
new_dg = lambda: mk(gE, nv, directed=True)
new_eu = lambda: mk(euE, ne)
new_c = lambda: nx.cycle_graph(range(1, 1001))
new_pair = lambda: (mk(gE, nv), mk(g2E, nv))
new_small_pair = lambda: (mk(hE, ns), mk(h2E, ns))

cold_bench("VertexDelete, 1000 of 100000 vertices (cold)", new_g, vdel)
check("VertexDelete, 1000 of 100000 vertices (cold)", cnt(vdel(new_g())))
cold_bench("EdgeAdd, 10000 edges (cold)", new_g, eadd)
check("EdgeAdd, 10000 edges (cold)", cnt(eadd(new_g())))
cold_bench("EdgeDelete, 10000 edges (cold)", new_g, edel)
check("EdgeDelete, 10000 edges (cold)", cnt(edel(new_g())))
cold_bench("Subgraph, 50000 of 100000 vertices (cold)", new_g, lambda h: h.subgraph(subV).copy())
check("Subgraph, 50000 of 100000 vertices (cold)", cnt(new_g().subgraph(subV).copy()))
cold_bench("GraphUnion, 2 x 100000 vertices (cold)", new_pair, lambda p: nx.compose(*p))
check("GraphUnion, 2 x 100000 vertices (cold)", cnt(nx.compose(*new_pair())))
cold_bench("GraphDifference, 2 x 2000 vertices (cold)", new_small_pair, lambda p: nx.difference(*p))
check("GraphDifference, 2 x 2000 vertices (cold)", cnt(nx.difference(*new_small_pair())))
cold_bench("GraphComplement, CycleGraph[1000] (cold)", new_c, nx.complement)
check("GraphComplement, CycleGraph[1000] (cold)", cnt(nx.complement(new_c())))
cold_bench("ReverseGraph, 199947 arcs (cold)", new_dg, lambda h: h.reverse(copy=True))
check("ReverseGraph, 199947 arcs (cold)", cnt(new_dg().reverse(copy=True)))
cold_bench("FindEulerianCycle, 100000 edges (cold)", new_eu, lambda h: list(nx.eulerian_circuit(h)))
check("FindEulerianCycle, 100000 edges (cold)", len(list(nx.eulerian_circuit(new_eu()))))
cold_bench("FindPath, 100000 vertices (cold)", new_g, lambda h: nx.shortest_path(h, 1, nv))
check("FindPath, 100000 vertices (cold)", 1)
