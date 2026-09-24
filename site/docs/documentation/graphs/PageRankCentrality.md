# PageRankCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PageRankCentrality[g, a] gives the PageRank of each vertex with damping factor a (default 0.85): x = a P^T x + (1-a)/n, where a vertex with no outgoing edge links to every vertex; the entries sum to 1. Edge weights are ignored.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= PageRankCentrality[CycleGraph[4]]
Out[1]= {0.25, 0.25, 0.25, 0.25}

In[2]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[2]= {0.213762, 0.264622, 0.307853, 0.213762}

In[3]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], 0.5]
Out[3]= {0.22449, 0.265306, 0.285714, 0.22449}

In[4]:= PageRankCentrality[StarGraph[5]]
Out[4]= {0.475676, 0.131081, 0.131081, 0.131081, 0.131081}

In[5]:= PageRankCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], 2]
Out[5]= PageRankCentrality[Graph[<4 vertices, 4 edges>], 2]
```

## Implementation notes

- Solves `x = a P^T x + (1 - a)/n`; dangling vertices jump uniformly;
  normalized to `Total[x] = 1`. Default `a = 0.85`; requires `0 <= a <= 1`
  (otherwise unevaluated). Weights ignored.
- Computed by (Jacobi) iteration with a correct stopping rule, converged to
  `1e-14` (Wolfram stops near `1e-9`, so they agree to ~9 digits).

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
