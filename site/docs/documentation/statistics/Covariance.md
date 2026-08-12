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

## Examples (3)

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
- For vectors, the unbiased estimate $\hat{\sigma}_{vw} = \frac{1}{n-1}\sum_i (v_i - \hat{\mu}_v)\overline{(w_i - \hat{\mu}_w)}$; the conjugate is on the **second** argument, so exact / complex / symbolic inputs yield exact / complex / symbolic output.
- For matrices, element $(i,j)$ is the covariance of column $i$ of `a` with column $j$ of `b`; `Covariance[a]` is symmetric.
- NDArray / packed real data uses a threaded centered inner product (vectors) or a BLAS gram (matrices); an integer sample degrades to the exact `List` path. Lowered inside `Compile[]`.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Covariance[]` reports `Covariance::argb`.

**Attributes:** `Protected`.

## See also

[List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
