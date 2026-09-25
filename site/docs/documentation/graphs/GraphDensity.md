# GraphDensity

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphDensity[g] gives the number of edges of g divided by the number of possible edges: (directed edges + 2 undirected edges)/(n(n-1)).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= GraphDensity[PathGraph[4]]
Out[1]= 1/2

In[2]:= GraphDensity[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= 1/3

In[3]:= GraphDensity[CompleteGraph[5]]
Out[3]= 1

In[4]:= GraphDensity[Graph[{1<->2, 2->3}]]
Out[4]= 1/2

In[5]:= GraphDensity[Graph[{1},{}]]
Out[5]= GraphDensity[Graph[<1 vertex, 0 edges>]]
```

## Implementation notes

- Computed as `(directed edges + 2 undirected edges)/(n (n - 1))`, exact.
- Weights ignored.
- Unevaluated for `n < 2`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
