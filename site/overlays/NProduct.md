### Worked examples

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5
```

```mathematica
In[1]:= NProduct[Cos[1/n], {n, 1, Infinity}]
Out[1]= 0.388536
```

### Notes

`NProduct[f, {i, imin, imax}]` numerically evaluates a product, with `imax`
allowed to be `Infinity`. The first case is a telescoping product:
`Product[1 - 1/n^2, {n, 2, Infinity}] = 1/2` exactly. The second converges to
about `0.388536`. `NProduct` is evaluated internally as `Exp[NSum[Log[f], ...]]`,
so the Euler–Maclaurin / Wynn's-epsilon machinery and convergence test of `NSum`
carry over. With `VerifyConvergence -> True` (default) a divergent product gives
`ComplexInfinity`. Use `WorkingPrecision` for arbitrary precision.
