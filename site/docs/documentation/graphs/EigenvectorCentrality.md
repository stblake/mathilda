# EigenvectorCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EigenvectorCentrality[g] gives the eigenvector centrality of each vertex, computed per strongly connected component (each component of k > 1 vertices weighted by k - 1, total 1; isolated vertices 0). EigenvectorCentrality[g, "In"] (default) uses incoming edges, [g, "Out"] outgoing ones. Edge weights are ignored.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]
Out[1]= {0.333333, 0.333333, 0.333333, 0.0}

In[2]:= EigenvectorCentrality[Graph[{1<->2,2<->3,3<->1,3<->4}]]
Out[2]= {0.269594, 0.269594, 0.315449, 0.145362}

In[3]:= EigenvectorCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], "Out"]
Out[3]= {0.333333, 0.333333, 0.333333, 0.0}

In[4]:= EigenvectorCentrality[Graph[{1->2,2->3}]]
Out[4]= {0.0, 0.0, 0.0}

In[5]:= EigenvectorCentrality[Graph[Join[EdgeList[CompleteGraph[4]], {5<->6,6<->7,5<->7}, {8<->9}]]]
Out[5]= {0.125, 0.125, 0.125, 0.125, 0.111111, 0.111111, 0.111111, 0.0833333, 0.0833333}
```

## Implementation notes

- Computed per strongly connected component: each component `C` with
  `|C| > 1` gets its Perron vector scaled to total
  `(|C| - 1)/Σ(|C'| - 1)`; single-vertex components get 0, so a DAG gives all
  zeros (*reverse-engineered*: reproduces Wolfram exactly on every
  disconnected / non-strongly-connected case tried, e.g. `K4 ⊔ K3 ⊔ K2` totals
  3:2:1).
- Weights ignored.
- **Algorithm:** restarted Arnoldi (Krylov dimension 40, BLAS
  re-orthogonalization, LAPACK on the Hessenberg matrix) per block,
  residual-tested to `1e-13`.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/gmet_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_init.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)
