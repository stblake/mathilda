# Correlation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Correlation[v, w]`**

gives the correlation between the vectors v and w, Covariance\[v, w\] / (StandardDeviation\[v\] StandardDeviation\[w\]).

**`Correlation[a, b]`**

gives the p\*q cross-correlation matrix between the columns of the matrices a and b.

**`Correlation[a]`**

gives the auto-correlation matrix of the columns of the matrix a; it is symmetric with a unit diagonal.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Correlation[{5, 3/4, 1}, {2, 1/2, 1}]
Out[1]= 2 Sqrt[3/13]

In[2]:= Correlation[{1.5, 3, 5, 10}, {2, 1.25, 15, 8}]
Out[2]= 0.475976

In[3]:= Correlation[{{a, b}, {c, d}}][[1, 1]]
Out[3]= 1
```

### Applications (4)

The correlation of two integer vectors is exact

```mathematica
In[4]:= Correlation[{1, 3, 5, 7}, {2, 3, 8, 9}]
Out[4]= 13/Sqrt[185]
```

A machine-real pair gives a machine-real correlation

```mathematica
In[5]:= Correlation[{1.5, 3., 5., 10.}, {2., 1.25, 15., 8.}]
Out[5]= 0.475976
```

A vector is perfectly correlated with itself

```mathematica
In[6]:= Correlation[{1, 3, 5, 7}, {1, 3, 5, 7}]
Out[6]= 1
```

One matrix gives the symmetric auto-correlation with a unit diagonal

```mathematica
In[7]:= Correlation[{{1, 2}, {3, 5}, {5, 7}}]
Out[7]= {{1, 5/2 Sqrt[3/19]}, {5/2 Sqrt[3/19], 1}}
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

**Algorithm.** `builtin_correlation` dispatches through `corrcov_dispatch(res, 1,
"Correlation")` — the same engine as `Covariance` with the `correlate` flag set —
accepting one to three arguments (`Correlation::argb` otherwise) and routing an
`NDArray` / packed argument to `nd_correlation` (`src/linalg/ndcorrcov.c`). A
correlation is a covariance normalized by the standard deviations:

- **`Correlation[v, w]`** — `g_corr_vectors` computes `Covariance[v, w]` and
  divides by `StandardDeviation[v] · StandardDeviation[w]` (the `n−1` factors
  cancel). For real data `−1 ≤ ρ ≤ 1`; exact / complex / symbolic data are carried
  through exactly as in `Covariance`.
- **`Correlation[a, b]`** — the `p×q` cross-correlation of the columns of two
  matrices; `g_matrix` divides each covariance cell by the two column standard
  deviations (precomputed per column in `sds_a` / `sds_b`).
- **`Correlation[a]`** — the `p×p` auto-correlation, symmetric with a unit
  diagonal. The diagonal is set outright rather than computed: it is the `Real`
  `1.` when the data has a machine real (`common_has_machine_real`), keeping the
  matrix one uniform type so `SymmetricMatrixQ` and `==` stay decidable, and the
  exact `Integer` `1` for exact / symbolic data.

**Data structures.** Shares `Covariance`'s `g_cov_vectors` / `g_matrix`; adds
`Expr**` arrays of per-column standard deviations for the normalization. Results
are `List`-of-`List` for the matrix forms.

**Complexity / limits.** `O(n)` for a vector pair; `O(p·q·n)` for the matrix forms,
plus `O(p + q)` standard-deviation reductions. On the `AWARE` list in `src/pack.c`
but not `INT64_OK` (shares `Covariance`'s `Rational`-degradation on an integer
buffer). Lowered inside `Compile[]`. Stays unevaluated for a single vector,
mismatched shapes, or fewer than two observations.

- `Protected`.
- A normalized covariance, $\rho_{vw} = \sigma_{vw} / (\sigma_v\,\sigma_w)$ with $\sigma_{vw} = \mathtt{Covariance}[v,w]$ and $\sigma_v = \mathtt{StandardDeviation}[v]$; $-1 \le \rho_{vw} \le 1$ for real data.
- The auto-correlation matrix `Correlation[a]` is symmetric with a unit diagonal (exact `1` for exact/symbolic data, `1.` for real data).
- Shares `Covariance`'s NDArray / packed / `Compile[]` fast paths.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Correlation[]` reports `Correlation::argb`.

**Attributes:** `Protected`.

## References

**See also:** [Covariance](../../statistics/Covariance/)

- Source: [`src/stats/corrcov.c`](https://github.com/stblake/mathilda/blob/main/src/stats/corrcov.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`Correlation` is a normalized `Covariance`:
`ρ = Covariance[v, w] / (StandardDeviation[v] StandardDeviation[w])`, with the
`n − 1` factors cancelling. For real data `−1 ≤ ρ ≤ 1`, and a vector is perfectly
(`ρ = 1`) correlated with itself.

`Correlation[a, b]` is the `p×q` cross-correlation of the columns of two matrices;
`Correlation[a]` is the `p×p` auto-correlation, symmetric with a unit diagonal. The
diagonal is the exact `1` for exact / symbolic data and the `Real` `1.` for real
data, so the matrix keeps one uniform type and `SymmetricMatrixQ` and `==` stay
decidable. It shares `Covariance`'s `NDArray` / packed and `Compile[]` fast paths,
stays unevaluated for a single vector, mismatched shapes, or fewer than two
observations, and `Correlation[]` reports `Correlation::argb`.
