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

### Notes

`Series[f, {x, x0, n}]` builds a power-series expansion to order `(x - x0)^n`, returning a `SeriesData` object that prints with a trailing `O[x]^(n+1)` term. It handles Taylor, Laurent (negative powers), and Puiseux (fractional powers) cases, as well as expansion around `Infinity` via the `x -> 1/u` substitution. Apply `Normal` to drop the order term and recover an ordinary polynomial, as in the `Exp` example above. `Series` is `HoldAll`, so the expansion variable is held unevaluated while the expansion point and order are read off.
