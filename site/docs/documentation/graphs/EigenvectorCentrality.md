# EigenvectorCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EigenvectorCentrality[g] gives the eigenvector centrality of each vertex, computed per strongly connected component (each component of k > 1 vertices weighted by k - 1, total 1; isolated vertices 0). EigenvectorCentrality[g, "In"] (default) uses incoming edges, [g, "Out"] outgoing ones. Edge weights are ignored.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {0.333333, 0.333333, 0.333333, 0.0}

In[2]:= EigenvectorCentrality[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[2]= {0.269594, 0.269594, 0.315449, 0.145362}

In[3]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Out"]
Out[3]= {0.333333, 0.333333, 0.333333, 0.0}

In[4]:= EigenvectorCentrality[Graph[{1->2,2->3}]]
Out[4]= {0.0, 0.0, 0.0}

In[5]:= EigenvectorCentrality[Graph[Join[EdgeList[CompleteGraph[4]], {5<->6,6<->7,5<->7}, {8<->9}]]]
Out[5]= {0.125, 0.125, 0.125, 0.125, 0.111111, 0.111111, 0.111111, 0.0833333, 0.0833333}
```

### Applications (6)

Every vertex of a cycle is equivalent

```mathematica
In[6]:= EigenvectorCentrality[CycleGraph[5]]
Out[6]= {0.2, 0.2, 0.2, 0.2, 0.2}
```

The hub is twice as central as each leaf

```mathematica
In[7]:= EigenvectorCentrality[StarGraph[5]]
Out[7]= {0.333333, 0.166667, 0.166667, 0.166667, 0.166667}
```

The middle vertex scores highest

```mathematica
In[8]:= EigenvectorCentrality[PathGraph[{1, 2, 3}]]
Out[8]= {0.292893, 0.414214, 0.292893}
```

Scores from outgoing arcs

```mathematica
In[9]:= EigenvectorCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1], DirectedEdge[1, 3]}], "Out"]
Out[9]= {0.43016, 0.245122, 0.324718}
```

Scores from incoming arcs, the default

```mathematica
In[10]:= EigenvectorCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1], DirectedEdge[1, 3]}], "In"]
Out[10]= {0.324718, 0.245122, 0.43016}
```

Normalised to total one on a connected graph

```mathematica
In[11]:= Total[EigenvectorCentrality[GridGraph[{3, 3}]]]
Out[11]= 1.0
```

## Implementation notes

**Algorithm.** `builtin_eigenvector_centrality` accepts `[g]`, `[g, "In"]` or `[g, "Out"]`. It splits the graph into strongly connected components (`gmet_scc`) and takes the Perron vector of each non-trivial component: `x_v` is proportional to the sum over arcs `u -> v` for `"In"` (the default), or over arcs `v -> u` for `"Out"`. Each component's vector is scaled to total `(|C| - 1) / sum(|C'| - 1)` over components, and single-vertex components get `0`. This weighting is reverse-engineered from Mathematica 15 output rather than documented. On a connected undirected graph it reduces to the usual normalization to sum 1. Perron vectors come from restarted Arnoldi (Krylov dimension up to 40, full re-orthogonalization, LAPACK `dgeev` on the small Hessenberg matrix) until `||M x - t x|| < 1e-13 |t|`; components of at most 64 vertices are solved densely.

**Data structures.** Two CSR views are built from the `Graph[List, List]` tree, one of out-arcs for the component search and one pull view (reversed for `"In"`). A block-local CSR is extracted per component. The result is a packed machine-real vector, cached by expression, and `EdgeWeight` is ignored.

**Complexity / limits.** Each Arnoldi step costs `O(m_C)` plus `O(k^2 n_C)` for orthogonalization. A component that fails to converge leaves the call unevaluated rather than answered inaccurately. Without LAPACK a shifted power iteration with the same residual test is used.

- Computed per strongly connected component: each component `C` with
  `|C| > 1` gets its Perron vector scaled to total
  `(|C| - 1)/Σ(|C'| - 1)`; single-vertex components get 0, so a DAG gives all
  zeros (*reverse-engineered*: reproduces Wolfram exactly on every
  disconnected / non-strongly-connected case tried, e.g. `K4 ⊔ K3 ⊔ K2` totals
  3:2:1).
- Weights ignored.
- **Algorithm:** restarted Arnoldi (Krylov dimension 40, BLAS
  re-orthogonalization, LAPACK on the Hessenberg matrix) per block,
  residual-tested to `1e-13`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_spectral.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_spectral.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

A vertex scores highly when it is linked to by other highly scoring vertices. The vector is the Perron eigenvector of the adjacency matrix, normalised to sum 1 on a connected graph. The optional `"In"` (default) or `"Out"` argument chooses whether incoming or outgoing arcs carry the score on a directed graph.

Disconnected and not strongly connected graphs are handled per strongly connected component, and vertices alone in a component score 0. The result is a list of machine reals.
