#!/usr/bin/env python3
"""Experiment 92 -- Graph structural predicates and ordering (networkx column).

Same ten cases as ``graph_predicates.m``, same graphs, same order.  networkx is
pure Python, so -- as in experiment 29 -- a Mathilda win here is a weaker claim
than a win against a compiled library.

Idioms are the ones a networkx user would reach for: ``v in g`` for VertexQ,
``g.has_edge`` for EdgeQ, ``density == 1`` for CompleteGraphQ (exact for a
simple graph), and ``is_forest`` for an undirected AcyclicGraphQ.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import networkx as nx

from harness import bench, check, require

require(["networkx", "networkx:is_directed_acyclic_graph", "networkx:topological_sort",
         "networkx:is_tree", "networkx:is_bipartite"])

nv = 100000
dag = nx.DiGraph()
dag.add_nodes_from(range(1, nv + 1))
dag.add_edges_from([(i, i + 1) for i in range(1, nv)]
                   + [(i, i + 2 + (7 * i) % 50) for i in range(1, nv - 51)])
tree = nx.Graph()
tree.add_nodes_from(range(1, nv + 1))
tree.add_edges_from((i, i // 2) for i in range(2, nv + 1))
bip = nx.Graph()
bip.add_nodes_from(range(1, nv + 1))
bip.add_edges_from([(i, i + 1) for i in range(1, nv)]
                   + [(i, i + 3) for i in range(1, nv - 2, 2)])
kn = nx.complete_graph(range(1, 401))


def wl(b):
    return "True" if b else "False"


def wl_list(xs):
    return "{" + ", ".join(xs) + "}"


def pair(f, a, b):
    return wl_list([wl(f(a)), wl(f(b))])


D = nx.DiGraph
G = nx.Graph

bench("AcyclicGraphQ, 100000-vertex DAG", lambda: nx.is_directed_acyclic_graph(dag))
check("AcyclicGraphQ, 100000-vertex DAG",
      pair(nx.is_directed_acyclic_graph, D([(1, 2), (2, 3), (1, 3)]), D([(1, 2), (2, 3), (3, 1)])))

bench("TopologicalSort, 100000-vertex DAG", lambda: list(nx.topological_sort(dag)))
check("TopologicalSort, 100000-vertex DAG",
      wl_list(str(v) for v in nx.topological_sort(D([(3, 1), (1, 2), (2, 4)]))))

bench("TreeGraphQ, 100000-vertex tree", lambda: nx.is_tree(tree))
check("TreeGraphQ, 100000-vertex tree",
      pair(nx.is_tree, G([(1, 2), (1, 3)]), G([(1, 2), (2, 3), (3, 1)])))

bench("AcyclicGraphQ, 100000-vertex tree", lambda: nx.is_forest(tree))
check("AcyclicGraphQ, 100000-vertex tree",
      pair(nx.is_forest, G([(1, 2), (3, 4)]), G([(1, 2), (2, 3), (3, 1)])))

bench("BipartiteGraphQ, 100000 vertices", lambda: nx.is_bipartite(bip))
check("BipartiteGraphQ, 100000 vertices",
      pair(nx.is_bipartite, G([(1, 2), (2, 3), (3, 4), (4, 1)]), G([(1, 2), (2, 3), (3, 1)])))

bench("CompleteGraphQ, K_400", lambda: nx.density(kn) == 1)
check("CompleteGraphQ, K_400",
      pair(lambda g: nx.density(g) == 1, G([(1, 2), (2, 3), (3, 1)]), G([(1, 2), (2, 3)])))

bench("UndirectedGraphQ, 100000-vertex tree", lambda: not tree.is_directed())
check("UndirectedGraphQ, 100000-vertex tree",
      pair(lambda g: not g.is_directed(), G([(1, 2)]), D([(1, 2)])))

bench("EmptyGraphQ, 100000-vertex tree", lambda: nx.is_empty(tree))
e2 = G(); e2.add_nodes_from([1, 2])
check("EmptyGraphQ, 100000-vertex tree", pair(nx.is_empty, e2, G([(1, 2)])))

bench("VertexQ, 100000-vertex DAG", lambda: nv in dag)
check("VertexQ, 100000-vertex DAG",
      wl_list([wl(2 in D([(1, 2)])), wl(3 in D([(1, 2)]))]))

bench("EdgeQ, 100000-vertex DAG", lambda: dag.has_edge(nv - 1, nv))
check("EdgeQ, 100000-vertex DAG",
      wl_list([wl(D([(1, 2)]).has_edge(1, 2)), wl(D([(1, 2)]).has_edge(2, 1))]))


# ---- COLD cases (see graph_predicates.m) ---------------------------------
# networkx keeps no per-graph answer cache, but the fresh-graph protocol is
# mirrored exactly so the three columns time the same thing.
import time


def cold_bench(label, mk, f, reps=5):
    ts = []
    for _ in range(reps):
        gg = mk()
        t0 = time.perf_counter()
        f(gg)
        ts.append(time.perf_counter() - t0)
    print("BENCH\t%s\t%.3f" % (label, 1000.0 * min(ts)), flush=True)


def fresh(template):
    g = template.__class__()
    g.add_nodes_from(template.nodes())
    g.add_edges_from(template.edges())
    return g


new_dag = lambda: fresh(dag)
new_tree = lambda: fresh(tree)
new_bip = lambda: fresh(bip)
new_kn = lambda: fresh(kn)

cold_bench("AcyclicGraphQ, 100000-vertex DAG (cold)", new_dag, nx.is_directed_acyclic_graph)
check("AcyclicGraphQ, 100000-vertex DAG (cold)", wl(nx.is_directed_acyclic_graph(new_dag())))
cold_bench("TopologicalSort, 100000-vertex DAG (cold)", new_dag, lambda g: list(nx.topological_sort(g)))
check("TopologicalSort, 100000-vertex DAG (cold)", len(list(nx.topological_sort(new_dag()))))
cold_bench("TreeGraphQ, 100000-vertex tree (cold)", new_tree, nx.is_tree)
check("TreeGraphQ, 100000-vertex tree (cold)", wl(nx.is_tree(new_tree())))
cold_bench("AcyclicGraphQ, 100000-vertex tree (cold)", new_tree, nx.is_forest)
check("AcyclicGraphQ, 100000-vertex tree (cold)", wl(nx.is_forest(new_tree())))
cold_bench("BipartiteGraphQ, 100000 vertices (cold)", new_bip, nx.is_bipartite)
check("BipartiteGraphQ, 100000 vertices (cold)", wl(nx.is_bipartite(new_bip())))
cold_bench("CompleteGraphQ, K_400 (cold)", new_kn, lambda g: nx.density(g) == 1)
check("CompleteGraphQ, K_400 (cold)", wl(nx.density(new_kn()) == 1))
cold_bench("UndirectedGraphQ, 100000-vertex tree (cold)", new_tree, lambda g: not g.is_directed())
check("UndirectedGraphQ, 100000-vertex tree (cold)", wl(not new_tree().is_directed()))
cold_bench("EmptyGraphQ, 100000-vertex tree (cold)", new_tree, nx.is_empty)
check("EmptyGraphQ, 100000-vertex tree (cold)", wl(nx.is_empty(new_tree())))
cold_bench("VertexQ, 100000-vertex DAG (cold)", new_dag, lambda g: nv in g)
check("VertexQ, 100000-vertex DAG (cold)", wl(nv in new_dag()))
cold_bench("EdgeQ, 100000-vertex DAG (cold)", new_dag, lambda g: g.has_edge(nv - 1, nv))
check("EdgeQ, 100000-vertex DAG (cold)", wl(new_dag().has_edge(nv - 1, nv)))
