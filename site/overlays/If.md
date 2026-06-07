---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= If[3 > 2, yes, no]
Out[1]= yes
```

```mathematica
In[1]:= If[x > 0, pos, neg]
Out[1]= If[x > 0, pos, neg]
```

```mathematica
In[1]:= If[1 == 2, a, b, c]
Out[1]= b
```

```mathematica
In[1]:= If[PrimeQ[7], prime, composite]
Out[1]= prime
```

### Notes

`If[cond, t, f]` evaluates `t` when `cond` is `True` and `f` when it is `False`.
`If` has the `HoldRest` attribute, so only the selected branch is evaluated — the
other is never run. When the condition is neither `True` nor `False` (e.g. the
symbolic `x > 0`), `If` returns unevaluated rather than guessing. The optional
fourth argument `u` is returned for that indeterminate case in the four-argument
form `If[cond, t, f, u]`.
