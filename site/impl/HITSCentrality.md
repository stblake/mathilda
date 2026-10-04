---
references:
  - "J. M. Kleinberg, *Authoritative sources in a hyperlinked environment*, J. ACM **46** (1999) 604-632."
source: src/graph/gmet_spectral.c
---
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
