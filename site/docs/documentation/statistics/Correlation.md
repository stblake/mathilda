# Correlation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Correlation[v, w]
    gives the correlation between the vectors v and w, Covariance[v, w] / (StandardDeviation[v] StandardDeviation[w]).
Correlation[a, b]
    gives the p*q cross-correlation matrix between the columns of the matrices a and b.
Correlation[a]
    gives the auto-correlation matrix of the columns of the matrix a; it is symmetric with a unit diagonal.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Correlation[{5, 3/4, 1}, {2, 1/2, 1}]
Out[1]= 2 Sqrt[3/13]

In[2]:= Correlation[{1.5, 3, 5, 10}, {2, 1.25, 15, 8}]
Out[2]= 0.475976

In[3]:= Correlation[{{a, b}, {c, d}}][[1, 1]]
Out[3]= 1
```

## Implementation notes

- `Protected`.
- A normalized covariance, $\rho_{vw} = \sigma_{vw} / (\sigma_v\,\sigma_w)$ with $\sigma_{vw} = \mathtt{Covariance}[v,w]$ and $\sigma_v = \mathtt{StandardDeviation}[v]$; $-1 \le \rho_{vw} \le 1$ for real data.
- The auto-correlation matrix `Correlation[a]` is symmetric with a unit diagonal (exact `1` for exact/symbolic data, `1.` for real data).
- Shares `Covariance`'s NDArray / packed / `Compile[]` fast paths.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Correlation[]` reports `Correlation::argb`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
