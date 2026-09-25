# EdgeConnectivity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeConnectivity[g] gives the minimum total weight of edges whose removal disconnects g (strongly, for directed graphs); EdgeConnectivity[g, s, t] gives the minimum s-t edge cut weight. EdgeWeight is used when present, else each edge counts 1.`**

## Examples (4)

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

## Implementation notes

- `Protected`; unevaluated on a non-graph.
- Weighted by `EdgeWeight`; strong connectivity for directed graphs.
- Integer weights give an Integer; Rational or Real ones give a Real; computed
  exactly in int64 by the flow/cut engine (Dinic; Nagamochi-Ibaraki for the
  undirected global value; `2(n-1)` bounded flows for directed graphs).

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/galg_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/galg_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)
