### Worked examples

```mathematica
In[1]:= SingularValueDecomposition[{{2, 0}, {0, 3}}]
Out[1]= {{{0, 1}, {1, 0}}, {{3, 0}, {0, 2}}, {{0, 1}, {1, 0}}}
```

The middle factor is the diagonal matrix of singular values (here `3` and `2`, in descending order); the outer factors permute the axes so the larger value comes first.

```mathematica
In[1]:= With[{r = SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]}, Chop[r[[1]] . r[[2]] . Transpose[r[[3]]]]]
Out[1]= {{1.0, 2.0}, {3.0, 4.0}}
```

The defining identity `m == u . sigma . ConjugateTranspose[v]` reconstructs the original matrix exactly (up to floating-point noise removed by `Chop`).

```mathematica
In[1]:= SingularValueDecomposition[N[{{1, 2}, {3, 4}}]]
Out[1]= {{{-0.404554, -0.914514}, {-0.914514, 0.404554}}, {{5.46499, 0.0}, {0.0, 0.365966}}, {{-0.576048, 0.817416}, {-0.817416, -0.576048}}}
```

### Notes

`SingularValueDecomposition[m]` returns `{u, sigma, v}` with `m == u . sigma . ConjugateTranspose[v]`, where `u` and `v` have orthonormal columns and `sigma` is diagonal with the singular values in descending order. Exact integer/rational matrices stay exact (singular values appear as `Sqrt[...]` forms when irrational); machine-precision real and complex inputs route through LAPACK, and high-precision MPFR input uses a one-sided Jacobi SVD. A two-argument form `SingularValueDecomposition[m, k]` keeps only the `k` largest singular values, and `SingularValueDecomposition[{m, a}]` computes the generalized SVD.
