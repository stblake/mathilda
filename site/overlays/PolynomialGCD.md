---
status: Stable
references:
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\" (3rd ed.), Ch. 6 & 11 (Euclidean and modular GCD)."
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (1992), Ch. 7 (polynomial GCD computation)."
---
### Worked examples

```mathematica
In[1]:= PolynomialGCD[x^2 - 1, x^2 + 2 x + 1]
Out[1]= 1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^4 - 1, x^2 - 1]
Out[1]= -1 + x^2
```

```mathematica
In[1]:= PolynomialGCD[x^2 - 1, x - 1]
Out[1]= -1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^3 - x, x^2 - x]
Out[1]= -x + x^2
```

```mathematica
In[1]:= PolynomialGCD[x^6 - 1, x^4 - 1, x^9 - 1]
Out[1]= -1 + x
```

```mathematica
In[1]:= PolynomialGCD[x^4 - 2, x^2 - Sqrt[2], Extension -> Sqrt[2]]
Out[1]= -Sqrt[2] + x^2
```

### Notes

`PolynomialGCD` returns the greatest common divisor of its polynomial
arguments, here over the rationals. The result is the highest-degree common
factor — for `x^2 - 1` and `(x+1)^2` that shared factor is `1 + x`, and for
`x^4 - 1` and `x^2 - 1` it is the full `x^2 - 1`. The output is normalized in
canonical term order and is not forced monic, so a shared `x` factor surfaces
as `-x + x^2`. The `Extension -> alpha` option computes the GCD over `Q(alpha)`
for `Sqrt[c]`, `c^(1/n)`, or `I`; the default treats any algebraic numbers as
independent variables.
