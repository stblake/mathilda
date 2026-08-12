#!/usr/bin/env python3
"""Experiment 29 -- Graph operations (networkx column).

Same seven kernels as ``graph_ops.m``, same order and sizes.

NOT A COMPILED BASELINE.  networkx is pure Python, so a Mathilda win here is a
weaker claim than a win against scipy or numpy.  What this row set is really for
is finding src/graph/ builtins that are ABSENT or accidentally quadratic --
nothing in the existing corpus has ever timed them (experiment 12 deliberately
avoided a graph type in order to isolate the gather).

The graph is deterministic -- a directed cycle plus 3i chords on 20000 vertices --
so all three systems build the same object and the checks are comparable.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import networkx as nx

from harness import bench, check, require

require(["networkx", "networkx:DiGraph", "networkx:shortest_path"])

nv = 20000
edges = ([(i, i % nv + 1) for i in range(1, nv + 1)]
         + [(i, (3 * i) % nv + 1) for i in range(1, nv + 1)])


def build():
    g = nx.DiGraph()
    g.add_edges_from(edges)
    return g


g = build()
tri = nx.DiGraph([(1, 2), (2, 3), (3, 1)])

bench("Graph construction, 40000 edges", build, reps=1)
check("Graph construction, 40000 edges", len(edges))

# Wolfram's VertexDegree is total degree (in + out) for a directed graph.
bench("VertexDegree, 20000 vertices", lambda: [d for _, d in g.degree()])
check("VertexDegree, 20000 vertices", sum(d for _, d in tri.degree()))

bench("EdgeCount", lambda: g.number_of_edges())
check("EdgeCount", tri.number_of_edges())

# Wolfram's ConnectedComponents on a directed graph gives WEAKLY connected ones.
bench("ConnectedComponents",
      lambda: list(nx.weakly_connected_components(g)), reps=1)
check("ConnectedComponents",
      nx.number_connected_components(nx.Graph([(1, 2), (3, 4)])))

bench("GraphDistance from vertex 1",
      lambda: nx.single_source_shortest_path_length(g, 1), reps=1)
check("GraphDistance from vertex 1",
      nx.shortest_path_length(nx.DiGraph([(1, 2), (2, 3)]), 1, 3))

bench("FindShortestPath 1 to 10000",
      lambda: nx.shortest_path(g, 1, 10000), reps=1)
check("FindShortestPath 1 to 10000",
      len(nx.shortest_path(nx.DiGraph([(1, 2), (2, 3)]), 1, 3)))

bench("AdjacencyMatrix, 20000 vertices",
      lambda: nx.adjacency_matrix(g), reps=1)
check("AdjacencyMatrix, 20000 vertices",
      int(nx.adjacency_matrix(tri).sum()))
