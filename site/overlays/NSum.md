### Worked examples

```mathematica
In[1]:= NSum[1/n^2, {n, 1, Infinity}]
Out[1]= 1.64493
```

```mathematica
In[1]:= NSum[(-1)^(n+1)/n, {n, 1, Infinity}]
Out[1]= 0.693147
```

```mathematica
In[1]:= NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[1]= 1.644934066848226436472415166649
```

```mathematica
In[1]:= NSum[1/n^4, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[1]= 1.082323233711138191516003696543
```

### Notes

`NSum[f, {i, imin, imax}]` numerically sums a series, with `imax` allowed to be
`Infinity`. The first two cases recover the Basel sum `Pi^2/6 = 1.64493...` and
the alternating harmonic sum `Log[2] = 0.693147...`. With `WorkingPrecision -> 30`
the Basel sum is computed to 30 digits, and `Sum[1/n^4]` returns
`Pi^4/90 = 1.082323233711...`. `Method -> Automatic` chooses Euler–Maclaurin for
monotone series, the Cohen–Villegas–Zagier method for alternating series, and
Wynn's epsilon otherwise. With `VerifyConvergence -> True` (default) a divergent
sum gives `ComplexInfinity`.
