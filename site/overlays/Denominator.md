---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational normal forms."
---
### Worked examples

```mathematica
In[1]:= Denominator[6/8]
Out[1]= 4
```

```mathematica
In[1]:= Denominator[(x+1)/(x-1)]
Out[1]= -1 + x
```

```mathematica
In[1]:= Denominator[a/b + c/d]
Out[1]= 1
```

```mathematica
In[1]:= Denominator[Together[a/b + c/d]]
Out[1]= b d
```

```mathematica
In[1]:= Denominator[(x^2-1)/((x-2)^3 (x+5))]
Out[1]= (5 + x) (-2 + x)^3
```

```mathematica
In[1]:= Denominator[x^(-2) y^3 z^(-1)]
Out[1]= x^2 z
```

### Notes

Denominator returns the bottom of the structural rational form. Rational constants
are reduced first, so `Denominator[6/8] = 4`. A non-fractional expression such as
an integer has denominator `1`. For symbolic quotients it returns the literal
denominator, e.g. `-1 + x` for `(x+1)/(x-1)`. Because Mathilda does not implicitly
`Together` a sum, `Denominator[a/b + c/d]` returns `1` (the expression is a `Plus`
with no overall denominator); call `Together` first to obtain `b d`.
