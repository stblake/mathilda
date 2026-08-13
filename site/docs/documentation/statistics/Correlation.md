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

## Examples (3)

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

- `Protected`.
- A normalized covariance, $\rho_{vw} = \sigma_{vw} / (\sigma_v\,\sigma_w)$ with $\sigma_{vw} = \mathtt{Covariance}[v,w]$ and $\sigma_v = \mathtt{StandardDeviation}[v]$; $-1 \le \rho_{vw} \le 1$ for real data.
- The auto-correlation matrix `Correlation[a]` is symmetric with a unit diagonal (exact `1` for exact/symbolic data, `1.` for real data).
- Shares `Covariance`'s NDArray / packed / `Compile[]` fast paths.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Correlation[]` reports `Correlation::argb`.

**Attributes:** `Protected`.

## References

**See also:** [Covariance](../../statistics/Covariance/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
