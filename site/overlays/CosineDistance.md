### Worked examples

```mathematica
In[1]:= CosineDistance[{1, 0}, {1, 1}]  (* 45 degrees apart: 1 - 1/Sqrt[2] *)
```

```mathematica
In[1]:= CosineDistance[{1, 2, 3}, {1, 2, 3}]  (* parallel vectors give 0 *)
```

```mathematica
In[1]:= CosineDistance[{0, 0}, {1, 2}]  (* a zero vector gives 0 by convention *)
```

### Notes

`CosineDistance[u, v]` is `1 - (u . Conjugate[v]) / (Norm[u] Norm[v])`, the
angular distance between two vectors. It runs over `[0, 2]`: `0` for parallel,
`1` for orthogonal, `2` for antiparallel. The `Conjugate` makes it correct for
complex vectors and is a no-op on reals.

Unlike the Euclidean family this is **not** a metric — it ignores magnitude and
violates the triangle inequality — so there is no squared form that ranks
identically. A zero vector on either side gives `0` (the quotient would
otherwise be the indeterminate `0/0`). Exact input stays exact and symbolic
vectors pass through.
