---
source: src/vectoranal.c
---
**Algorithm.** `builtin_curl` (through the `vecop` front end) computes the
Cartesian curl as a generalized Levi-Civita contraction:
`(1/k!) Sum eps_{a..ij..} d_{x_i} f_{j..}`, where `k = ncube_depth(f)` is the
field's depth and the result has depth `n - k - 1`. `curl_perm_recur` enumerates
all `n!` permutations; `curl_emit_perm` splits each as `(a_part | i | j_part)`,
reads the leaf `f_{j_part}`, and accumulates `sign * D[f_{j_part}, x_i]`
(`perm_sign` is `(-1)^inversions`) into the cell indexed by `a_part`.
`curl_build_nested` then folds each cell's terms into `Plus[...]/k!` and
assembles the depth-`(n-k-1)` nested list, reduced with one `eval_and_free`.
This covers the three familiar cases uniformly: a 2-D vector gives a scalar, a
3-D vector gives a vector, and a rank-2 tensor gives a scalar.

**Data structures.** A `curl_ctx` holds the field, the variables, and a
per-cell array of growable term buffers (`curl_cell_append` doubles capacity);
`Expr` trees are built with `mk_d`/`mk_neg`/`mk_fnN_adopt`. There is no ND or
`Compile[]` path — the output is a symbolic derivative. The 3-argument
curvilinear form (`curl_chart`, `n = 2` or `3` only) builds the orthonormal-basis
curl from the chart's Lamé factors, e.g. the 3-D component
`(1/(h_j h_k))[D[h_k f_k, x_j] - D[h_j f_j, x_k]]` cyclically.

**Complexity / limits.** The permutation enumeration is `O(n!)`, so the
Cartesian path is bounded to `2 <= n <= 6`; outside that, or when the field is a
scalar (`k < 1`) or `k > n - 1` (negative result depth), it returns `NULL`. The
chart form declines (`NULL`) for `n` other than 2 or 3, a non-vector field, or
an unrecognised chart (emitting `Curl::chart`) — tensor curl in a curvilinear
basis needs a metric and is out of scope.
