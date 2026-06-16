### Worked examples

```mathematica
In[1]:= N[GoldenRatio]
Out[1]= 1.61803
```

```mathematica
In[1]:= N[GoldenRatio, 40]
Out[1]= 1.6180339887498948482045868343656381177203
```

```mathematica
In[1]:= N[GoldenRatio^2 - GoldenRatio - 1, 40]
Out[1]= 0.0
```

```mathematica
In[1]:= FromContinuedFraction[{1, {1}}]
Out[1]= 1/2 (1 + Sqrt[5])
```

```mathematica
In[1]:= Round[N[(GoldenRatio^15 - (1 - GoldenRatio)^15)/Sqrt[5]]]
Out[1]= 610
```

### Notes

`GoldenRatio` is `phi = (1 + Sqrt[5])/2`, the positive root of `x^2 == x + 1`;
evaluating that polynomial at `phi` numerically returns `0.0`, confirming the
defining relation. It has the simplest possible continued fraction `[1; 1, 1, ...]`,
so `FromContinuedFraction[{1, {1}}]` recovers it exactly. Through Binet's formula
`Fibonacci[n] == (phi^n - (1 - phi)^n)/Sqrt[5]`, the closed form at `n = 15`
rounds to `610 = Fibonacci[15]`. It is a protected `Constant` (so
`D[GoldenRatio, x]` is `0`) evaluated to any precision by `N`.
