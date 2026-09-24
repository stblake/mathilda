# GlobalClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GlobalClusteringCoefficient[g] gives 3 x (number of triangles) / (number of connected triples) of g, exactly.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= GlobalClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 3/5

In[2]:= GlobalClusteringCoefficient[PathGraph[4]]
Out[2]= 0

In[3]:= GlobalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 3/4

In[4]:= GlobalClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= GlobalClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

## Implementation notes

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- Undirected: `3T/Σ C(d, 2)`, with `T` the number of triangles and the sum
  over the vertex degrees `d`.
- Directed (*reverse-engineered*): triangles are directed 3-cycles, with the
  same in/out/reciprocal counting as `LocalClusteringCoefficient`.

**Attributes:** `Protected`.

## References

**See also:** [LocalClusteringCoefficient](../../graphs/LocalClusteringCoefficient/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
