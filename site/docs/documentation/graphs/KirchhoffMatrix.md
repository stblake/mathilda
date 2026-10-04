# KirchhoffMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KirchhoffMatrix[g] gives the Kirchhoff (Laplacian) matrix D - A of g, with D the diagonal matrix of vertex degrees (incident edges) and A the adjacency matrix. Dense (Wolfram returns a SparseArray); weights are ignored.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= KirchhoffMatrix[PathGraph[3]]
Out[1]= {{1, -1, 0}, {-1, 2, -1}, {0, -1, 1}}

In[2]:= KirchhoffMatrix[Graph[{1->2,2->3}]]
Out[2]= {{1, -1, 0}, {0, 2, -1}, {0, 0, 1}}

In[3]:= NDArrayQ[KirchhoffMatrix[PathGraph[3]]]
Out[3]= True

In[4]:= KirchhoffMatrix[5]
Out[4]= KirchhoffMatrix[5]
```

### Options (1)

```mathematica
In[5]:= KirchhoffMatrix[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[5]= {{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}}
```

### Applications (5)

Degree 2 on the diagonal, -1 for each neighbour

```mathematica
In[6]:= KirchhoffMatrix[CycleGraph[4]]
Out[6]= {{2, -1, 0, -1}, {-1, 2, -1, 0}, {0, -1, 2, -1}, {-1, 0, -1, 2}}
```

The hub has degree 3, the leaves degree 1

```mathematica
In[7]:= KirchhoffMatrix[StarGraph[4]]
Out[7]= {{3, -1, -1, -1}, {-1, 1, 0, 0}, {-1, 0, 1, 0}, {-1, 0, 0, 1}}
```

A directed edge fills only one off-diagonal entry

```mathematica
In[8]:= KirchhoffMatrix[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3]}]]
Out[8]= {{1, -1, 0}, {0, 2, -1}, {0, 0, 1}}
```

Undirected rows sum to zero

```mathematica
In[9]:= Total[KirchhoffMatrix[CompleteGraph[5]], 2]
Out[9]= 0
```

The Laplacian is singular

```mathematica
In[10]:= Det[KirchhoffMatrix[CycleGraph[4]]]
Out[10]= 0
```

## Implementation notes

**Algorithm.** `builtin_kirchhoff_matrix` forms the Laplacian `D - A` in a single pass over the edge list. For each edge `k` with endpoints `a = eu[k]`, `b = ev[k]` it increments both diagonal entries, so each vertex's diagonal entry is the number of edges incident to it. It then sets `K[a][b] = -1`, and also `K[b][a] = -1` when the edge is undirected. Edge weights are ignored. A directed edge `u -> v` therefore contributes to both degrees but only to the `u`-to-`v` off-diagonal entry.

**Data structures.** The edge-index view from `graph_edge_indices` supplies the parallel arrays `eu`, `ev` and `directed`. The result is a dense `n*n` `int64_t` buffer returned as a packed int64 `NDArray` (`ndbuild_open`), with a plain nested `List` of integers as the fallback. Mathematica returns a `SparseArray`; Mathilda has none, so this is its `Normal` form.

**Complexity / limits.** `O(n^2)` time and space for the zero-fill, plus `O(m)` for the scan. Repeated edges accumulate on the diagonal but not off it. The empty graph and non-graph arguments leave the call unevaluated.

- `D` is the diagonal matrix of the number of incident edges of each vertex,
  `A` the (directed) adjacency matrix; weights ignored.
- **Deviation:** returns a dense packed Integer matrix (Wolfram returns a
  `SparseArray`; this is its `Normal`).

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [SparseArray](../../data-structures/SparseArray/), [Normal](../../data-structures/Normal/)

- Source: [`src/graph/gmet_distance.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_distance.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The result is `D - A`, where `D` holds the number of edges incident to each vertex and `A` is the adjacency matrix. Edge weights are ignored.

Mathematica returns a `SparseArray`; Mathilda returns the equivalent dense integer matrix, which is what `Normal` of Mathematica's answer gives.
