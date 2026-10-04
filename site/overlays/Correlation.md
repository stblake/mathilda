### Worked examples

```mathematica
In[1]:= Correlation[{1, 3, 5, 7}, {2, 3, 8, 9}]  (* the correlation of two integer vectors is exact *)
Out[1]= 13/Sqrt[185]
```

```mathematica
In[1]:= Correlation[{1.5, 3., 5., 10.}, {2., 1.25, 15., 8.}]  (* a machine-real pair gives a machine-real correlation *)
Out[1]= 0.475976
```

```mathematica
In[1]:= Correlation[{1, 3, 5, 7}, {1, 3, 5, 7}]  (* a vector is perfectly correlated with itself *)
Out[1]= 1
```

```mathematica
In[1]:= Correlation[{{1, 2}, {3, 5}, {5, 7}}]  (* one matrix gives the symmetric auto-correlation with a unit diagonal *)
Out[1]= {{1, 5/2 Sqrt[3/19]}, {5/2 Sqrt[3/19], 1}}
```

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
