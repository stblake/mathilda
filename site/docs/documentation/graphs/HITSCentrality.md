# HITSCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HITSCentrality[g] gives {h, a}: h is the principal eigenvector of A^T A (per block of vertices sharing in-neighbours, each block of k > 1 vertices weighted k - 1, total 1) and a = A h, with A the adjacency matrix (Wolfram's convention). Edge weights are ignored.`**

## Examples (5)

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

### Applications (1)

{hubs, authorities}

```mathematica
In[5]:= HITSCentrality[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]
Out[5]= {{0.0, 0.381966, 0.618034, 0.0}, {1.0, 0.618034, 0.0, 0.0}}
```

## Implementation notes

**Algorithm.** `builtin_hits_centrality` returns `{h, a}`, the hub and authority vectors of
Kleinberg's HITS. Rather than the classic alternating power iteration on `A` and `A^T`, it
computes the hub vector as the principal (Perron) eigenvector of `A^T A` and then sets the
authority vector `a = A h`. `A^T A` is block-diagonal over classes of vertices that *share an
in-neighbour* — found by union-find over the co-citation relation — so each block is solved
independently and the blocks are recombined. Each block's Perron vector is found densely via
LAPACK `dgeev` when it is small (≤ 64 vertices), and by a restarted Arnoldi iteration (Krylov
dimension up to 40, residual tolerance `1e-13`) for larger blocks, falling back to shifted
power iteration where LAPACK is unavailable. Each block of `k` vertices is scaled to total mass
`(k-1)/Σ(size-1)`, so singleton (trivial) blocks contribute 0; this reproduces Wolfram's
otherwise-undocumented normalisation.

**Data structures.** Its own out-arc `GmetCSR`, a union-find `par[]` with class map `cls[]`, the
iteration vectors `h` and `a`, a local `BlockOp` CSR per block, and, on the Arnoldi path, a
Krylov basis `V` with an upper-Hessenberg `H`. `EdgeWeight` is ignored.

**Complexity / limits.** Dominated by the per-block eigensolve — `O(k^3)` dense, or ≤ 40
matrix-vector products per Arnoldi restart (≤ 2000 restarts). The empty graph returns `{}`; a
block that fails to converge returns unevaluated; a mixed or non-graph argument returns
unevaluated.

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

- J. M. Kleinberg, *Authoritative sources in a hyperlinked environment*, J. ACM **46** (1999) 604-632.
- Source: [`src/graph/gmet_spectral.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_spectral.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The result is the pair `{h, a}`: `h` the hub scores and `a` the authority scores, each a vector
in `VertexList` order. A vertex is a good hub when it points to good authorities, and a good
authority when it is pointed to by good hubs.

`EdgeWeight` is ignored — HITS is a purely structural measure. The empty graph returns `{}`, and
a graph whose vertices never share an in-neighbour (so every co-citation class is a singleton)
scores all zeros.
