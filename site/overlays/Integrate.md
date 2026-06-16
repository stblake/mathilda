---
status: Partial
references:
  - "Bronstein, \"Symbolic Integration I: Transcendental Functions\", 2nd ed. (Springer, 2005)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 11–12."
---
### Worked examples

```mathematica
In[1]:= Integrate[1/(1 + x^2), x]
Out[1]= ArcTan[x]
```

```mathematica
In[1]:= Integrate[1/x, x]
Out[1]= Log[x]
```

```mathematica
In[1]:= Integrate[Cos[x], x]
Out[1]= Sin[x]
```

```mathematica
In[1]:= Integrate[x^3 + x, x]
Out[1]= 1/2 x^2 + 1/4 x^4
```

```mathematica
In[1]:= Integrate[1/(x^3 + 1), x]
Out[1]= 1/3 Log[1 + x] + ArcTan[(-1 + 2 x)/Sqrt[3]]/Sqrt[3] - 1/6 Log[1 - x + x^2]
```

```mathematica
In[1]:= Integrate[(x^2 + 1)/(x^4 + 1), x]
Out[1]= ArcTan[x/Sqrt[2]]/Sqrt[2] + ArcTan[(x + x^3)/Sqrt[2]]/Sqrt[2]
```

```mathematica
In[1]:= Integrate[x*Exp[x], x]
Out[1]= -E^x + x E^x

In[2]:= Integrate[1/(x*Log[x]), x]
Out[2]= Log[Log[x]]
```

### Notes

`Integrate[f, x]` computes the indefinite integral via a cascade: Bronstein's rational-function algorithm, then the Risch–Norman (`pmint`) heuristic, then the lazy-loaded CRC integral tables; `Method -> "<name>"` pins a single subroutine. Antiderivatives are returned without an integration constant and are not always simplified — for example `Integrate[Sin[x], x]` returns `-(1 + Cos[x])` rather than `-Cos[x]`. Definite integration is **not** supported in this build: `Integrate[x^2, {x, 0, 1}]` threads `Integrate` over the bound list and returns the garbage form `{1/3 x^3, Integrate[x^2, 0], Integrate[x^2, 1]}` instead of `1/3`. Restrict use to the indefinite, single-variable form.
