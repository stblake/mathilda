---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", on polynomial GCD computation."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), on rational function simplification."
---
### Worked examples

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 + 2 x + 1)/(x + 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x^2 - 2 x + 1)]
Out[1]= (1 + x)/(-1 + x)
```

```mathematica
In[1]:= Cancel[(x^4 - 1)/(x^2 - 1)]
Out[1]= 1 + x^2
```

```mathematica
In[1]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[1]= Sqrt[2] + x
```

### Notes

Cancel removes the common polynomial factors between numerator and denominator by
dividing out their GCD, without performing a partial-fraction split. It reduces
`(x^2-1)/(x-1)` to `1 + x` and recognises perfect-square numerators such as
`(x^2+2x+1)/(x+1)`. When a common factor remains after cancellation, the result
stays a reduced quotient: `(x^2-1)/(x^2-2x+1)` becomes `(1+x)/(-1+x)` because both
share the factor `(x-1)` but the leftover `(x+1)/(x-1)` is already in lowest
terms. Cancel does not expand the surviving factors back out.

With the `Extension -> alpha` option, cancellation is performed over the
algebraic field `Q(alpha)`: `(x^2 - 2)/(x - Sqrt[2])` factors as
`(x - Sqrt[2])(x + Sqrt[2])` over `Q(Sqrt[2])`, so the denominator divides out
and the quotient collapses to `Sqrt[2] + x`. The default `Extension -> None`
treats algebraic numbers as opaque and would leave that quotient intact.
