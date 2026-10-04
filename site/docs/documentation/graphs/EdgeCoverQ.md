# EdgeCoverQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeCoverQ[g, elist] gives True if elist is a set of edges of g touching every vertex of g.`**

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[4,3]}]
Out[1]= True

In[2]:= IndependentEdgeSetQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[2]= False

In[3]:= IndependentEdgeSetQ[Graph[{1->2,3->4}], {DirectedEdge[2,1]}]
Out[3]= False

In[4]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[3,4]}]
Out[4]= True

In[5]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,2], UndirectedEdge[2,3]}]
Out[5]= False

In[6]:= EdgeCoverQ[CycleGraph[4], {UndirectedEdge[1,3]}]
Out[6]= False
```

### Applications (6)

Two disjoint edges touch all four vertices

```mathematica
In[7]:= EdgeCoverQ[CycleGraph[4], {1 <-> 2, 3 <-> 4}]
Out[7]= True
```

Vertices 3 and 4 are missed

```mathematica
In[8]:= EdgeCoverQ[CycleGraph[4], {1 <-> 2}]
Out[8]= False
```

An undirected edge matches in either orientation

```mathematica
In[9]:= EdgeCoverQ[PathGraph[{1, 2, 3}], {2 <-> 1, 2 <-> 3}]
Out[9]= True
```

One spoke leaves the others uncovered

```mathematica
In[10]:= EdgeCoverQ[StarGraph[4], {UndirectedEdge[1, 2]}]
Out[10]= False
```

A pair that is not an edge of the graph gives false

```mathematica
In[11]:= EdgeCoverQ[CycleGraph[4], {1 <-> 3}]
Out[11]= False
```

A non-graph gives false

```mathematica
In[12]:= EdgeCoverQ[x, {}]
Out[12]= False
```

## Implementation notes

**Algorithm.** `EdgeCoverQ[g, es]` is a linear check, not a search (`gm_edge_set_q` with `cover = 1`). It gives `False` unless `g` is a valid graph and `es` is a list. Each element must be an edge of `g`: `gm_edge_of` accepts `UndirectedEdge`/`TwoWayRule` (matching in either orientation) and `DirectedEdge`/`Rule` (only as given), resolved through `graph_has_edge`. Every endpoint is marked in a visited mask, and the answer is `True` exactly when every vertex of `g` was touched. Its sibling `IndependentEdgeSetQ` shares the code and instead rejects an edge whose end was already marked.

**Data structures.** A byte array `hit[n]` indexed by vertex position, with positions looked up by `galg_vertex_arg` from the vertex list of the `Graph[List, List]` tree.

**Complexity / limits.** `O(n + |es|)` plus the edge-membership lookups. Never stays unevaluated for two arguments: a non-graph, a non-list, or an element that is not an edge of `g` gives `False`. A graph with an isolated vertex has no edge cover. Repeated edges are harmless.

- `Protected` membership predicates; give `False` when `g` is not a graph.
- An element that is not an edge of `g` gives `False`. An `UndirectedEdge`
  matches either orientation; a `DirectedEdge` matches only as given.

**Attributes:** `Protected`.

## References

**See also:** [IndependentEdgeSetQ](../../graphs/IndependentEdgeSetQ/)

- Source: [`src/graph/galg_mis.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_mis.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

The second argument is a list of edges of the graph, written `a <-> b` or `UndirectedEdge[a, b]` for undirected edges and `a -> b` for directed ones; a directed edge only matches as given. The predicate is `True` exactly when every vertex of the graph is an endpoint of one of the listed edges.

It is a linear check, never a search: use it to validate a candidate cover. A list element that is not an edge of the graph, or a non-list, gives `False`.
