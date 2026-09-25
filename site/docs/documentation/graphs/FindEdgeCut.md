# FindEdgeCut

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindEdgeCut[g] gives a minimum set of edges whose removal disconnects g (strongly, for directed graphs); FindEdgeCut[g, s, t] gives a minimum s-t edge cut. Uses EdgeWeight as capacity when present. Edges are returned in EdgeList order.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindEdgeCut[CycleGraph[5]]
Out[1]= {1 <-> 2, 5 <-> 1}

In[2]:= FindEdgeCut[CycleGraph[6], 1, 4]
Out[2]= {1 <-> 2, 6 <-> 1}

In[3]:= FindEdgeCut[Graph[{1->2,2->3,3->1}]]
Out[3]= {1 -> 2}

In[4]:= FindEdgeCut[x]
Out[4]= FindEdgeCut[x]
```

### Options (1)

```mathematica
In[5]:= FindEdgeCut[Graph[{1,2,3},{UndirectedEdge[1,2],UndirectedEdge[2,3],UndirectedEdge[3,1]}, EdgeWeight->{2,3,4}]]
Out[5]= {1 <-> 2, 2 <-> 3}
```

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Cuts are weighted by `EdgeWeight`; the edges are returned in `EdgeList` order.
- The `s`-`t` cut is the one closest to `s`, as in Mathematica.
- Directed graphs use strong connectivity.
- Uses the same exact int64 flow machinery as `FindMaximumFlow` /
  `FindMinimumCut` (Dinic; Nagamochi-Ibaraki for undirected global cuts).

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/), [EdgeList](../../graphs/EdgeList/), [FindMaximumFlow](../../graphs/FindMaximumFlow/), [FindMinimumCut](../../graphs/FindMinimumCut/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
