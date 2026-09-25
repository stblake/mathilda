# KirchhoffMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KirchhoffMatrix[g] gives the Kirchhoff (Laplacian) matrix D - A of g, with D the diagonal matrix of vertex degrees (incident edges) and A the adjacency matrix. Dense (Wolfram returns a SparseArray); weights are ignored.`**

## Examples (5)

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

## Implementation notes

- `D` is the diagonal matrix of the number of incident edges of each vertex,
  `A` the (directed) adjacency matrix; weights ignored.
- **Deviation:** returns a dense packed Integer matrix (Wolfram returns a
  `SparseArray`; this is its `Normal`).

**Attributes:** `Protected`.

## References

**See also:** [D](../../calculus/D/), [Normal](../../data-structures/Normal/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
