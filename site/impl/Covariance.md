---
source: src/stats/corrcov.c
---
**Algorithm.** `builtin_covariance` dispatches through `corrcov_dispatch(res, 0,
"Covariance")`, which accepts one to three arguments (`Covariance::argb`
otherwise) and routes an `NDArray` / packed argument to the buffer fast path
`nd_covariance` (`src/linalg/ndcorrcov.c`). The three shapes:

- **`Covariance[v, w]`** — a scalar, when both are `VectorQ`. Handled by
  `g_cov_vectors`: both must be equal-length `List`s with `n ≥ 2`. A machine-double
  fast path (all elements real-numeric, at least one inexact, none complex) computes
  the two means and returns `(1/(n−1)) Σ (vᵢ − μ_v)(wᵢ − μ_w)` as a `Real`.
  Otherwise the unbiased estimate `(1/(n−1)) Σ (vᵢ − μ_v) Conjugate[wᵢ − μ_w]` is
  built as an expression and evaluated — so exact input stays exact, complex stays
  complex, symbolic stays symbolic, with no `int64` overflow risk. The `Conjugate`
  is on the **second** argument.
- **`Covariance[a, b]`** — the `p×q` cross-covariance of the columns of two
  `n`-row matrices (`g_matrix`). It transposes each operand, requires `≥ 2`
  observations, and fills cell `(i, j)` with `g_cov_vectors(col_i(a), col_j(b))`; a
  degenerate column yields `0` to keep the shape rectangular.
- **`Covariance[a]`** — the `p×p` auto-covariance `Covariance[a, a]`, via the same
  `g_matrix` with the `auto_form` flag.

**Data structures.** `Expr` trees driven through `eval_and_free`; the vector path
uses an `Expr**` of `n` product terms before the `Plus`. The matrix path keeps the
two transposes as `List`-of-`List` and assembles an `Expr**` of rows. The buffer
path is a threaded centered inner product off the packed `double` buffer (vectors)
or a BLAS gram `A_cᵀ B_c` (matrices).

**Complexity / limits.** `O(n)` for a vector pair; `O(p·q·n)` for the matrix forms
(a BLAS `O(p q n)` gram on the buffer). On the `AWARE` list in `src/pack.c` but not
`INT64_OK` — an integer sample's covariance is a `Rational` no float64 slot holds,
so an integer buffer degrades to the exact `List` path, like `Variance`. Lowered
inside `Compile[]`. Stays unevaluated for a single vector, mismatched shapes, or
fewer than two observations.
