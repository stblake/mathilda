# GraphCenter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphCenter[g] gives the vertices of g with minimum eccentricity; {} unless g is (strongly) connected.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GraphCenter[PathGraph[5]]
Out[1]= {3}

In[2]:= GraphPeriphery[PathGraph[5]]
Out[2]= {1, 5}

In[3]:= GraphCenter[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= {}

In[4]:= GraphPeriphery[Graph[{1<->2,3<->4}]]
Out[4]= {}
```

### Options (1)

```mathematica
In[5]:= GraphCenter[Graph[{1,2,3},{1->2,2->3}, EdgeWeight->{1,2}]]
Out[5]= {1}
```

## Implementation notes

- *(w)* weight-aware.
- Unweighted and not strongly connected (not connected, if undirected): `{}`.
- Weighted: taken over the weighted eccentricities, so a weighted digraph that
  is not strongly connected can have a non-empty center
  (*reverse-engineered*).
- Reduced from the cached MS-BFS per-source summary (see `GraphDistanceMatrix`).

**Attributes:** `Protected`.

## References

**See also:** [GraphPeriphery](../../graphs/GraphPeriphery/), [GraphDistanceMatrix](../../graphs/GraphDistanceMatrix/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
