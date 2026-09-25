# MeanGraphDistance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanGraphDistance[g] gives the mean distance over all ordered pairs of distinct vertices of g: exact for unweighted graphs, a machine real for weighted ones, Infinity unless g is (strongly) connected.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= MeanGraphDistance[PathGraph[4]]
Out[1]= 5/3

In[2]:= MeanGraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= Infinity

In[3]:= MeanGraphDistance[Graph[{1},{}]]
Out[3]= MeanGraphDistance[Graph[<1 vertex, 0 edges>]]
```

### Options (1)

```mathematica
In[4]:= MeanGraphDistance[Graph[{1,2,3},{1<->2,2<->3,1<->3}, EdgeWeight->{1,1,5}]]
Out[4]= 1.33333
```

## Implementation notes

- *(w)* weight-aware (machine real when weighted).
- Averages over ordered pairs of distinct vertices; exact when unweighted.
- `Infinity` when some pair is unreachable (in particular, unweighted and not
  strongly connected).
- Unevaluated for a single vertex.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
