# MeanClusteringCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanClusteringCoefficient[g] gives the mean of the local clustering coefficients of g, exactly.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= MeanClusteringCoefficient[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[1]= 7/12

In[2]:= MeanClusteringCoefficient[CompleteGraph[4]]
Out[2]= 1

In[3]:= MeanClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[3]= 5/8

In[4]:= MeanClusteringCoefficient[Graph[{1<->2, 2->3}]]
Out[4]= MeanClusteringCoefficient[Graph[<3 vertices, 2 edges>]]
```

## Implementation notes

- Exact, weights ignored, mixed graphs unevaluated (as in Wolfram).
- The arithmetic mean of `LocalClusteringCoefficient[g]` (undirected or
  directed definition accordingly).

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
