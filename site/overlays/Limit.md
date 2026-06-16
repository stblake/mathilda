---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 3."
---
### Worked examples

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1
```

```mathematica
In[1]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[1]= 2
```

```mathematica
In[1]:= Limit[(1 + a/x)^x, x -> Infinity]
Out[1]= E^a
```

```mathematica
In[1]:= Limit[(Sin[x] - x + x^3/6)/x^5, x -> 0]
Out[1]= 1/120
```

```mathematica
In[1]:= Limit[(x^x - x)/(1 - x + Log[x]), x -> 1]
Out[1]= -2
```

```mathematica
In[1]:= Limit[x - Sqrt[x^2 + x], x -> Infinity]
Out[1]= -1/2
```

```mathematica
In[1]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[1]= 5
```

```mathematica
In[1]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[1]= Infinity

In[2]:= Limit[1/x, x -> 0, Direction -> "FromBelow"]
Out[2]= -Infinity
```

### Notes

`Limit[f, x -> a]` resolves the standard removable-singularity and indeterminate forms, including the classic `(1 + 1/x)^x -> E` and `0/0` cancellations such as `(x^2 - 1)/(x - 1)`. The `Direction` option selects one-sided (`"FromAbove"`/`"FromBelow"`) or complex approaches; the default is two-sided. Results may be a finite value, `Infinity`, `ComplexInfinity`, `Indeterminate`, an `Interval`, or the original expression unevaluated when the limit cannot be determined. Iterated and joint multivariate limits are supported through the list forms of the second argument.
