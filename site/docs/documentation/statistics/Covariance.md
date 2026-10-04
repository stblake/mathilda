# Covariance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Covariance[v, w]`**

gives the unbiased covariance estimate between the vectors v and w, (1/(n-1)) Sum\[(v\_i - Mean\[v\]) Conjugate\[w\_i - Mean\[w\]\]\].

**`Covariance[a, b]`**

gives the p\*q cross-covariance matrix between the columns of the matrices a and b.

**`Covariance[a]`**

gives the auto-covariance matrix of the columns of the matrix a, i.e. Covariance\[a, a\].

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Covariance[{1, 3/2}, {2, 11}]
Out[1]= 9/4

In[2]:= Covariance[{2 + I, 3 - 2 I, 5 + 4 I}, {I, 1 + 2 I, 10 - 5 I}]
Out[2]= -7/3 + 56/3*I

In[3]:= Covariance[{{1, 2}, {3, 4}, {5, 7}}]
Out[3]= {{4, 5}, {5, 19/3}}
```

### Applications (5)

The unbiased covariance of two integer vectors is exact

```mathematica
In[4]:= Covariance[{1, 3, 5, 7}, {2, 3, 8, 9}]
Out[4]= 26/3
```

A machine-real pair gives a machine-real covariance

```mathematica
In[5]:= Covariance[{1.5, 3., 5., 10.}, {2., 1.25, 15., 8.}]
Out[5]= 11.2604
```

A vector with itself is its variance

```mathematica
In[6]:= Covariance[{1, 3, 5, 7}, {1, 3, 5, 7}]
Out[6]= 20/3
```

Complex data: the conjugate falls on the second argument

```mathematica
In[7]:= Covariance[{2 + I, 3 - 2 I, 5 + 4 I}, {I, 1 + 2 I, 10 - 5 I}]
Out[7]= -7/3 + 56/3*I
```

One matrix gives the auto-covariance of its columns

```mathematica
In[8]:= Covariance[{{1, 2}, {3, 4}, {5, 7}}]
Out[8]= {{4, 5}, {5, 19/3}}
```

## Algorithm

corrcov.c -- Covariance[] and Correlation[].

```text
Covariance[v, w]   covariance between two length-n vectors (a scalar)
Covariance[a, b]   p x q cross-covariance of the columns of two n-row matrices
Covariance[a]      p x p auto-covariance of a matrix, i.e. Covariance[a, a]
Correlation[...]   the same three shapes, normalized by the standard deviations
```

For length-n vectors the covariance is

```text
  (1/(n-1)) Sum_i (v_i - Mean[v]) Conjugate[w_i - Mean[w]]
```

(the conjugate is on the SECOND argument), and the correlation divides that by StandardDeviation[v] StandardDeviation[w] (the (n-1) factors cancel). The matrix forms apply the vector definition to each pair of columns.

Following variance.c, the exact/complex/symbolic work is built as sub-expressions and evaluated, so exact input yields exact output, complex yields complex, and symbolic yields symbolic — with no int64-overflow risk. A fast machine-double path covers real numeric vectors; an NDArray / packed-array argument takes the buffer fast path in src/linalg/ndcorrcov.c.

See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

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

- `Protected`.
- For vectors, the unbiased estimate $\hat{\sigma}_{vw} = \frac{1}{n-1}\sum_i (v_i - \hat{\mu}_v)\overline{(w_i - \hat{\mu}_w)}$; the conjugate is on the **second** argument, so exact / complex / symbolic inputs yield exact / complex / symbolic output.
- For matrices, element $(i,j)$ is the covariance of column $i$ of `a` with column $j$ of `b`; `Covariance[a]` is symmetric.
- NDArray / packed real data uses a threaded centered inner product (vectors) or a BLAS gram (matrices); an integer sample degrades to the exact `List` path. Lowered inside `Compile[]`.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Covariance[]` reports `Covariance::argb`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/)

- Source: [`src/stats/corrcov.c`](https://github.com/stblake/mathilda/blob/main/src/stats/corrcov.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

For two vectors `Covariance[v, w]` is the unbiased estimate
`(1/(n − 1)) Σ (vᵢ − Mean[v]) Conjugate[wᵢ − Mean[w]]`. The conjugate is on the
**second** argument, so `Covariance[v, v]` reproduces `Variance[v]` exactly, and a
complex pair returns a complex covariance. Exact input stays exact and symbolic
stays symbolic — the expression is built and evaluated rather than reduced to
machine doubles, so there is no `int64`-overflow risk.

`Covariance[a, b]` is the `p×q` cross-covariance of the columns of two `n`-row
matrices, and `Covariance[a]` the `p×p` auto-covariance `Covariance[a, a]`. A real
`NDArray` / packed buffer uses a threaded centered inner product (vectors) or a
BLAS gram (matrices); an integer sample degrades to the exact `List` path. The
head stays unevaluated for a single vector, mismatched shapes, or fewer than two
observations; `Covariance[]` reports `Covariance::argb`.
