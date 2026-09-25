# HITSCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HITSCentrality[g] gives {h, a}: h is the principal eigenvector of A^T A (per block of vertices sharing in-neighbours, each block of k > 1 vertices weighted k - 1, total 1) and a = A h, with A the adjacency matrix (Wolfram's convention). Edge weights are ignored.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= HITSCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {{0.5, 0.0, 0.0, 0.5}, {0.0, 0.0, 1.0, 0.0}}

In[2]:= HITSCentrality[Graph[{1->2,1->3,2->3}]]
Out[2]= {{0.0, 0.381966, 0.618034}, {1.0, 0.618034, 0.0}}

In[3]:= HITSCentrality[PathGraph[5]]
Out[3]= {{0.166667, 0.166667, 0.333333, 0.166667, 0.166667}, {0.166667, 0.5, 0.333333, 0.5, 0.166667}}

In[4]:= HITSCentrality[Graph[{1<->2, 3->4}]]
Out[4]= {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}
```

## Implementation notes

- Returns `{h, A.h}`; `h` is built like `EigenvectorCentrality` but for
  `A^T A`, whose blocks are the classes of vertices sharing an in-neighbour
  (*reverse-engineered*; reproduces Wolfram exactly, including degenerate
  spectra such as `PathGraph[5]` and the all-zero answer for mixed graphs whose
  classes are singletons).
- Weights ignored. Spectral blocks use the same restarted Arnoldi solver as
  `EigenvectorCentrality`.

**Attributes:** `Protected`.

## References

**See also:** [EigenvectorCentrality](../../graphs/EigenvectorCentrality/)

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
