---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= Apply[Plus, {1, 2, 3, 4}]
Out[1]= 10
```

```mathematica
In[1]:= f @@ {a, b, c}
Out[1]= f[a, b, c]
```

```mathematica
In[1]:= Apply[List, a + b + c]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Apply[f, {{a, b}, {c, d}}, {1}]
Out[1]= {f[a, b], f[c, d]}
```

```mathematica
In[1]:= Apply[Times, Range[10]]
Out[1]= 3628800
```

```mathematica
In[1]:= Apply[GCD, {84, 126, 210}]
Out[1]= 42

In[2]:= Apply[Plus, Table[1/k^2, {k, 1, 6}]]
Out[2]= 5369/3600
```

### Notes

`f @@ expr` is the shorthand for `Apply[f, expr]`: it replaces the head of `expr`
with `f`. `Apply[Plus, list]` is the standard idiom for summing a list, since
the list's `List` head is swapped for `Plus`. Because addition is stored as
`Plus[...]`, `Apply[List, a + b + c]` recovers the summands. A level
specification like `{1}` applies the head replacement to each element at that
level instead of the whole expression. Folding a variadic head over a generated
list is a common pattern: `Apply[Times, Range[10]]` is `10!`, `Apply[GCD, ...]`
collapses a list to a single greatest common divisor, and `Apply[Plus, ...]`
over exact rationals returns the exact partial sum `5369/3600`.
