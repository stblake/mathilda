### Worked examples

```mathematica
In[1]:= SquaredEuclideanDistance[{0, 0}, {3, 4}]  (* the squared 3-4-5 distance *)
```

```mathematica
In[1]:= SquaredEuclideanDistance[{1, 2, 2}, {0, 0, 0}]  (* 1 + 4 + 4 = 9 *)
```

```mathematica
In[1]:= SquaredEuclideanDistance[3, 8]  (* scalars act as 1-vectors *)
```

### Notes

`SquaredEuclideanDistance[u, v]` is `Sum Abs[u_i - v_i]^2` — the Euclidean
distance without the final square root. Because no root is taken, the result is
exact for exact input (`SquaredEuclideanDistance[{1/3, 0}, {0, 1/7}]` is
`58/441`, not a float), and since squaring is monotone on non-negatives, ranking
on the square orders points identically to ranking on the true distance. That is
exactly what lets `FindClusters` partition exact multi-dimensional data without
ever introducing an irrational.

Both arguments must be scalars or equal-length lists; complex components use
their modulus and symbolic input survives.
