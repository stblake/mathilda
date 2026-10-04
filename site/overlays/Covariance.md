### Worked examples

```mathematica
In[1]:= Covariance[{1, 3, 5, 7}, {2, 3, 8, 9}]  (* the unbiased covariance of two integer vectors is exact *)
Out[1]= 26/3
```

```mathematica
In[1]:= Covariance[{1.5, 3., 5., 10.}, {2., 1.25, 15., 8.}]  (* a machine-real pair gives a machine-real covariance *)
Out[1]= 11.2604
```

```mathematica
In[1]:= Covariance[{1, 3, 5, 7}, {1, 3, 5, 7}]  (* a vector with itself is its variance *)
Out[1]= 20/3
```

```mathematica
In[1]:= Covariance[{2 + I, 3 - 2 I, 5 + 4 I}, {I, 1 + 2 I, 10 - 5 I}]  (* complex data: the conjugate falls on the second argument *)
Out[1]= -7/3 + 56/3*I
```

```mathematica
In[1]:= Covariance[{{1, 2}, {3, 4}, {5, 7}}]  (* one matrix gives the auto-covariance of its columns *)
Out[1]= {{4, 5}, {5, 19/3}}
```

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
