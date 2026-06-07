---
source: src/list.c
---
**Algorithm.** `builtin_transpose` swaps the levels of a rectangular nested-`List` array. It
measures the array shape with `get_array_dimensions` (requiring depth ≥ 2 and rectangularity),
then either uses the default permutation `{2, 1, 3, …}` (one-argument form swaps the first two
levels) or the explicit permutation given as the second argument. `build_transposed` recursively
materialises the output array by mapping each output index path back to an input index path
through the permutation and copying the leaf via `get_element_at`. For a 2-D matrix (list of
rows) this is the ordinary `m[i][j] -> m[j][i]` swap. Returns `NULL` (unevaluated) for
non-rectangular or non-`List` inputs. `ConjugateTranspose` is `Conjugate[Transpose[...]]`.
