# Covariance

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Covariance[v, w]
    gives the unbiased covariance estimate between the vectors v and w, (1/(n-1)) Sum[(v_i - Mean[v]) Conjugate[w_i - Mean[w]]].
Covariance[a, b]
    gives the p*q cross-covariance matrix between the columns of the matrices a and b.
Covariance[a]
    gives the auto-covariance matrix of the columns of the matrix a, i.e. Covariance[a, a].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Covariance[{1, 3/2}, {2, 11}]
Out[1]= 9/4

In[2]:= Covariance[{2 + I, 3 - 2 I, 5 + 4 I}, {I, 1 + 2 I, 10 - 5 I}]
Out[2]= -7/3 + 56/3*I

In[3]:= Covariance[{{1, 2}, {3, 4}, {5, 7}}]
Out[3]= {{4, 5}, {5, 19/3}}
```

## Implementation notes

- `Protected`.
- For vectors, the unbiased estimate $\hat{\sigma}_{vw} = \frac{1}{n-1}\sum_i (v_i - \hat{\mu}_v)\overline{(w_i - \hat{\mu}_w)}$; the conjugate is on the **second** argument, so exact / complex / symbolic inputs yield exact / complex / symbolic output.
- For matrices, element $(i,j)$ is the covariance of column $i$ of `a` with column $j$ of `b`; `Covariance[a]` is symmetric.
- NDArray / packed real data uses a threaded centered inner product (vectors) or a BLAS gram (matrices); an integer sample degrades to the exact `List` path. Lowered inside `Compile[]`.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Covariance[]` reports `Covariance::argb`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
