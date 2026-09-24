# EccentricityCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EccentricityCentrality[g] gives 1/VertexEccentricity[g, v] for each vertex v (0 when the eccentricity is 0). Uses EdgeWeight as lengths.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= EccentricityCentrality[PathGraph[5]]
Out[1]= {0.25, 0.333333, 0.5, 0.333333, 0.25}

In[2]:= EccentricityCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.333333, 0.5, 0.5, 0.0}

In[3]:= EccentricityCentrality[Graph[{1,2},{}]]
Out[3]= {0.0, 0.0}
```

### Options (1)

```mathematica
In[4]:= EccentricityCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= {0.333333, 0.5, 0.0}
```

## Implementation notes

- *(w)* weight-aware; machine reals.
- `e(v)` is measured over the vertices reachable from `v` (also when
  weighted); the centrality is 0 when `e(v) = 0`.
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
