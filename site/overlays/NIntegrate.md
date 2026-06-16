### Worked examples

```mathematica
In[1]:= NIntegrate[Sin[x], {x, 0, Pi}]
Out[1]= 2.0
```

```mathematica
In[1]:= NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}]
Out[1]= 1.77245
```

```mathematica
In[1]:= NIntegrate[Exp[-x^2], {x, -Infinity, Infinity}, WorkingPrecision -> 30]
Out[1]= 1.772453850905516027298167483341
```

```mathematica
In[1]:= NIntegrate[Sin[x]/x, {x, 0, Infinity}]
Out[1]= 1.5708
```

```mathematica
In[1]:= NIntegrate[Log[x] Log[1 - x], {x, 0, 1}]
Out[1]= 0.355066
```

### Notes

`NIntegrate[f, {x, a, b}]` approximates a definite integral. The Gaussian
example reproduces `Sqrt[Pi] = 1.77245...`, computed to 30 digits with
`WorkingPrecision -> 30` via the double-exponential rule on the infinite range.
The Dirichlet integral `Sin[x]/x` over `[0, Infinity]` is the oscillatory case,
returning `Pi/2`. The final integral has the closed form `2 - Pi^2/6 =
0.355066...`. `Method -> Automatic` selects globally-adaptive Gauss-Kronrod,
double-exponential, Levin oscillatory, or Monte-Carlo schemes per region.
Endpoints may be infinite or complex (a contour).
