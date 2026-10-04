---
references:
  - "P. Erdős and A. Rényi, *On random graphs I*, Publ. Math. Debrecen **6** (1959) 290-297."
source: src/graph/generators.c
---
**Algorithm.** `builtin_random_graph` samples from the Erdős–Rényi `G(n, m)` model: a simple
undirected graph on vertices `1..n` with exactly `m` edges, every `m`-edge graph equally likely.
The `m` edges are drawn without replacement from the `n(n-1)/2` candidate pairs, but the
candidate list is never materialised: `random_sample_indices(maxe, m)` produces exactly the
draws `RandomSample` would (so it honours `SeedRandom`), and each sampled index is decoded to a
vertex pair by binary search over the triangular row offsets. The assembled `Graph[...]` is
re-validated by the evaluator's `builtin_graph`.

**Data structures.** Only the `m` sampled indices and the resulting edge `List`; the decode is
`O(m log n)` time and `O(m)` memory — no `O(n^2)` candidate array.

**Complexity / limits.** `O(m log n)`. Guards: `n` or `m` negative returns unevaluated; `n`
above `2^31-1` returns unevaluated (overflow); `m` greater than `n(n-1)/2` returns unevaluated.
`RandomGraph[{n, m}, k]` returns a list of `k` independent graphs (`k = 0` gives `{}`); on any
sub-failure the whole list is discarded. Always undirected.
