---
source: src/graph/gmet_spectral.c
---
**Algorithm.** `builtin_eigenvector_centrality` accepts `[g]`, `[g, "In"]` or `[g, "Out"]`. It splits the graph into strongly connected components (`gmet_scc`) and takes the Perron vector of each non-trivial component: `x_v` is proportional to the sum over arcs `u -> v` for `"In"` (the default), or over arcs `v -> u` for `"Out"`. Each component's vector is scaled to total `(|C| - 1) / sum(|C'| - 1)` over components, and single-vertex components get `0`. This weighting is reverse-engineered from Mathematica 15 output rather than documented. On a connected undirected graph it reduces to the usual normalization to sum 1. Perron vectors come from restarted Arnoldi (Krylov dimension up to 40, full re-orthogonalization, LAPACK `dgeev` on the small Hessenberg matrix) until `||M x - t x|| < 1e-13 |t|`; components of at most 64 vertices are solved densely.

**Data structures.** Two CSR views are built from the `Graph[List, List]` tree, one of out-arcs for the component search and one pull view (reversed for `"In"`). A block-local CSR is extracted per component. The result is a packed machine-real vector, cached by expression, and `EdgeWeight` is ignored.

**Complexity / limits.** Each Arnoldi step costs `O(m_C)` plus `O(k^2 n_C)` for orthogonalization. A component that fails to converge leaves the call unevaluated rather than answered inaccurately. Without LAPACK a shifted power iteration with the same residual test is used.
