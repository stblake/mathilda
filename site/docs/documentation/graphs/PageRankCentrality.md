# PageRankCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PageRankCentrality[g, a] gives the PageRank of each vertex with damping factor a (default 0.85): x = a P^T x + (1-a)/n, where a vertex with no outgoing edge links to every vertex; the entries sum to 1. Edge weights are ignored.`**

## Examples (10)

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

### Applications (5)

A symmetric graph gives equal ranks

```mathematica
In[6]:= PageRankCentrality[CycleGraph[4]]
Out[6]= {0.25, 0.25, 0.25, 0.25}
```

The hub collects the most rank at the default damping 0.85

```mathematica
In[7]:= PageRankCentrality[StarGraph[4]]
Out[7]= {0.47973, 0.173423, 0.173423, 0.173423}
```

Lower damping flattens the ranks toward uniform

```mathematica
In[8]:= PageRankCentrality[StarGraph[4], 1/2]
Out[8]= {0.416667, 0.194444, 0.194444, 0.194444}
```

The sink at the end of a directed path ranks highest

```mathematica
In[9]:= PageRankCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3]}]]
Out[9]= {0.184417, 0.341171, 0.474412}
```

Ranks always sum to one

```mathematica
In[10]:= Total[PageRankCentrality[PathGraph[{1, 2, 3, 4, 5}]]]
Out[10]= 1.0
```

## Implementation notes

**Algorithm.** `builtin_pagerank_centrality` accepts `[g]` or `[g, a]` with damping `0 <= a <= 1`, default `0.85`. It solves `x = a P^T x + (1 - a)/n` with `P` the row-stochastic walk matrix, by power iteration from the uniform vector. A vertex with no out-arc (dangling) jumps uniformly: each step adds `a * (sum of dangling mass) / n` to the teleport term. After every pull the vector is renormalized to sum 1 before the L1 change is measured, and iteration stops when that change falls below `1e-14` (or after 10000 steps). An undirected edge counts as an arc each way, and `EdgeWeight` is ignored.

**Data structures.** A CSR of in-arcs is built from the `Graph[List, List]` tree. Out-degrees are stored as reciprocals. Two `double` buffers are swapped each step, with a scratch vector holding `x[u] / outdeg[u]`, so each step is a pure gather over in-neighbours, split over the thread team for large graphs. The result is a packed machine-real vector, cached by expression.

**Complexity / limits.** `O(m + n)` per iteration, with the iteration count set by `a` (roughly `log(1e-14) / log(a)`). It declines for `a` outside `[0, 1]` or a non-numeric `a`.

- Solves `x = a P^T x + (1 - a)/n`; dangling vertices jump uniformly;
  normalized to `Total[x] = 1`. Default `a = 0.85`; requires `0 <= a <= 1`
  (otherwise unevaluated). Weights ignored.
- Computed by (Jacobi) iteration with a correct stopping rule, converged to
  `1e-14` (Wolfram stops near `1e-9`, so they agree to ~9 digits).

**Attributes:** `Protected`.

## References

- S. Brin and L. Page, *The anatomy of a large-scale hypertextual Web search engine*, Computer Networks and ISDN Systems **30** (1998) 107-117.
- Source: [`src/graph/gmet_spectral.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_spectral.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The optional second argument is the damping factor `a` with `0 <= a <= 1`; the default is `0.85`. The ranks solve `x = a P^T x + (1 - a)/n` and sum to 1. A vertex with no outgoing edge redistributes its rank uniformly over all vertices.

Iteration runs to an L1 change below `10^-14`, tighter than Mathematica's stopping rule, so the two agree to about nine digits. Edge weights are ignored.
