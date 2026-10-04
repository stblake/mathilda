---
references:
  - "S. Brin and L. Page, *The anatomy of a large-scale hypertextual Web search engine*, Computer Networks and ISDN Systems **30** (1998) 107-117."
source: src/graph/gmet_spectral.c
---
**Algorithm.** `builtin_pagerank_centrality` accepts `[g]` or `[g, a]` with damping `0 <= a <= 1`, default `0.85`. It solves `x = a P^T x + (1 - a)/n` with `P` the row-stochastic walk matrix, by power iteration from the uniform vector. A vertex with no out-arc (dangling) jumps uniformly: each step adds `a * (sum of dangling mass) / n` to the teleport term. After every pull the vector is renormalized to sum 1 before the L1 change is measured, and iteration stops when that change falls below `1e-14` (or after 10000 steps). An undirected edge counts as an arc each way, and `EdgeWeight` is ignored.

**Data structures.** A CSR of in-arcs is built from the `Graph[List, List]` tree. Out-degrees are stored as reciprocals. Two `double` buffers are swapped each step, with a scratch vector holding `x[u] / outdeg[u]`, so each step is a pure gather over in-neighbours, split over the thread team for large graphs. The result is a packed machine-real vector, cached by expression.

**Complexity / limits.** `O(m + n)` per iteration, with the iteration count set by `a` (roughly `log(1e-14) / log(a)`). It declines for `a` outside `[0, 1]` or a non-numeric `a`.
