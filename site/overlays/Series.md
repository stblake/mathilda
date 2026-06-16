---
status: Stable
references:
  - "Geddes, Czapor & Labahn, \"Algorithms for Computer Algebra\" (Kluwer, 1992), ch. 3."
---
### Worked examples

```mathematica
In[1]:= Series[Sin[x], {x, 0, 5}]
Out[1]= x - 1/6 x^3 + 1/120 x^5 + O[x]^6
```

```mathematica
In[1]:= Series[1/(1 - x), {x, 0, 4}]
Out[1]= 1 + x + x^2 + x^3 + x^4 + O[x]^5
```

```mathematica
In[1]:= Series[Log[1 + x], {x, 0, 4}]
Out[1]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + O[x]^5
```

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 3}]]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3
```

The expansion of `Tan` recovers the tangent numbers as coefficients:

```mathematica
In[1]:= Series[Tan[x], {x, 0, 7}]
Out[1]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7 + O[x]^8
```

A Laurent expansion appears automatically when `f` has a pole at the base point;
`1/(E^x - 1)` starts at `x^(-1)` and its coefficients are the Bernoulli numbers
over factorials:

```mathematica
In[1]:= Series[1/(Exp[x] - 1), {x, 0, 4}]
Out[1]= 1/x - 1/2 + 1/12 x - 1/720 x^3 + O[x]^5
```

`Series` also handles symbolic-exponent expansions such as `x^x`, where the
coefficients involve `Log[x]`:

```mathematica
In[1]:= Series[x^x, {x, 0, 3}]
Out[1]= 1 + Log[x] x + 1/2 Log[x]^2 x^2 + 1/6 Log[x]^3 x^3 + O[x]^4
```

Expansion at infinity recovers the classic limit `(1 + 1/x)^x -> E`, with the
correction terms made explicit:

```mathematica
In[1]:= Series[(1 + 1/x)^x, {x, Infinity, 2}]
Out[1]= E + (-1/2 E)/x + 11/24 E (1/x)^2 + O[1/x]^3
```

### Notes

`Series[f, {x, x0, n}]` builds a power-series expansion to order `(x - x0)^n`, returning a `SeriesData` object that prints with a trailing `O[x]^(n+1)` term. It handles Taylor, Laurent (negative powers), and Puiseux (fractional powers) cases, as well as expansion around `Infinity` via the `x -> 1/u` substitution. Apply `Normal` to drop the order term and recover an ordinary polynomial, as in the `Exp` example above. `Series` is `HoldAll`, so the expansion variable is held unevaluated while the expansion point and order are read off.
