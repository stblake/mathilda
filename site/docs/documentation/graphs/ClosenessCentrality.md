# ClosenessCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ClosenessCentrality[g] gives, for each vertex v, r/s where r is the number of vertices reachable from v and s the sum of their distances from v (0 if v reaches none). Uses EdgeWeight as lengths.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ClosenessCentrality[PathGraph[4]]
Out[1]= {0.5, 0.75, 0.75, 0.5}

In[2]:= ClosenessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.5, 0.6, 0.75, 0.0}

In[3]:= ClosenessCentrality[Graph[{1,2,3},{}]]
Out[3]= {0.0, 0.0, 0.0}
```

### Options (1)

```mathematica
In[4]:= ClosenessCentrality[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[4]= {0.5, 0.5, 0.0}
```

## Implementation notes

- *(w)* weight-aware; machine reals, packed.
- For a vertex `v`: `r/s`, with `r` the number of vertices reachable from `v`
  and `s` the sum of their distances; 0 if none are reachable.
- Reduced from the cached MS-BFS per-source summary (weighted: Dijkstra per
  source); see `GraphDistanceMatrix`.

**Attributes:** `Protected`.

## References

**See also:** [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
