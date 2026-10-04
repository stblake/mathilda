# EdgeConnectivity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeConnectivity[g] gives the minimum total weight of edges whose removal disconnects g (strongly, for directed graphs); EdgeConnectivity[g, s, t] gives the minimum s-t edge cut weight. EdgeWeight is used when present, else each edge counts 1.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EdgeConnectivity[CompleteGraph[5]]
Out[1]= 4

In[2]:= EdgeConnectivity[CycleGraph[6], 1, 4]
Out[2]= 2

In[3]:= EdgeConnectivity[Graph[{1->2,2->3,3->1}]]
Out[3]= 1
```

### Options (1)

```mathematica
In[4]:= EdgeConnectivity[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[3,1]}, EdgeWeight->{2,3,4}]]
Out[4]= 5
```

### Applications (4)

Two edges must go to break a cycle

```mathematica
In[5]:= EdgeConnectivity[CycleGraph[5]]
Out[5]= 2
```

Equals the minimum degree, 3

```mathematica
In[6]:= EdgeConnectivity[CompleteGraph[4]]
Out[6]= 3
```

A single bridge disconnects a path

```mathematica
In[7]:= EdgeConnectivity[PathGraph[{1, 2, 3}]]
Out[7]= 1
```

The minimum 1-2 edge cut

```mathematica
In[8]:= EdgeConnectivity[CompleteGraph[4], 1, 2]
Out[8]= 3
```

## Implementation notes

**Algorithm.** `builtin_edge_connectivity` gives the minimum total capacity of a
set of edges whose removal disconnects `g` (`EdgeWeight` as the capacity when
present, else `1` per edge). `EdgeConnectivity[g]` is the *global* minimum cut
and `EdgeConnectivity[g, s, t]` the minimum `s-t` cut. For an undirected graph
the global cut uses **Nagamochi-Ibaraki**: each round computes a
maximum-adjacency order, lowers the bound `lambda` to the least weighted degree,
and contracts every edge whose MA-attachment reaches `lambda` (so it cannot lie
on a lighter cut) — it returns Stoer-Wagner's value while usually contracting
most of the graph per round. A directed global cut is the minimum over `v != v0`
of the flows `v0 -> v` and `v -> v0`. The `s-t` form is a single max flow.

**Data structures.** Capacities are scaled to exact `int64` in a `GfCap` (integer
capacities stay exact; reals scale by a common power of two; `Infinity` gets a
finite stand-in above every finite sum). Flows run Dinic on a CSR residual
network (`GfNet`: `first[]`/`to[]`/`rev[]`/`res[]` arcs, reverse arcs adjacent),
with level-BFS and current-arc blocking-flow search. The undirected cut keeps its
own CSR plus a lazy max-heap for the MA order.

**Complexity / limits.** All arithmetic is on `int64`; the result is an `Integer`
for integer capacities and a `Real` otherwise. A graph with `< 2` vertices is
left unevaluated, as is a negative or symbolic capacity.

- `Protected`; unevaluated on a non-graph.
- Weighted by `EdgeWeight`; strong connectivity for directed graphs.
- Integer weights give an Integer; Rational or Real ones give a Real; computed
  exactly in int64 by the flow/cut engine (Dinic; Nagamochi-Ibaraki for the
  undirected global value; `2(n-1)` bounded flows for directed graphs).

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/)

- H. Nagamochi and T. Ibaraki, *Computing edge-connectivity in multigraphs and capacitated graphs*, SIAM J. Discrete Math. **5** (1992) 54-66.
- E. A. Dinic, *Algorithm for solution of a problem of maximum flow in a network with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280.
- Source: [`src/graph/galg_flow.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_flow.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`EdgeConnectivity[g]` is the weight of a global minimum edge cut — the cheapest
set of edges whose removal disconnects the graph (strongly, for a directed
graph). `EdgeConnectivity[g, s, t]` restricts this to cuts that separate `s`
from `t`, which by the max-flow/min-cut theorem equals the maximum `s-t` flow
under unit (or `EdgeWeight`) capacities.

With integer weights the answer is an exact integer. A graph with fewer than
two vertices is left unevaluated.
