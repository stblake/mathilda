### Worked examples

```mathematica
In[1]:= HypergeometricPFQ[{1, 1}, {2}, z]
Out[1]= -Log[1 - z]/z
```

Common upper and lower parameters cancel, and any of the specialised forms
(`0F1`, `1F1`, `2F1`) is just a particular shape of `HypergeometricPFQ`. A
higher `3F2` evaluates numerically by direct summation:

```mathematica
In[1]:= N[HypergeometricPFQ[{1, 2, 3}, {4, 5}, 1/2], 40]
Out[1]= 1.1898747542564229318256831180919799547257
```

A non-positive integer upper parameter terminates the series to a polynomial:

```mathematica
In[1]:= HypergeometricPFQ[{-3, 1}, {1}, z]
Out[1]= 1 - 3 z + 3 z^2 - z^3
```

### Notes

`HypergeometricPFQ[{a1, ...}, {b1, ...}, z]` is the generalized hypergeometric function `pFq`, the series `Sum[(Product[Pochhammer[a_i, k]] / Product[Pochhammer[b_j, k]]) z^k / k!]`. It converges for all `z` when `p <= q`, and for `|z| < 1` when `p == q + 1`; a non-positive integer upper parameter truncates it to a polynomial. The specialised heads `Hypergeometric0F1`, `Hypergeometric1F1`, and `Hypergeometric2F1` are convenience wrappers.
