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

### Notes

`f @@ expr` is the shorthand for `Apply[f, expr]`: it replaces the head of `expr`
with `f`. `Apply[Plus, list]` is the standard idiom for summing a list, since
the list's `List` head is swapped for `Plus`. Because addition is stored as
`Plus[...]`, `Apply[List, a + b + c]` recovers the summands. A level
specification like `{1}` applies the head replacement to each element at that
level instead of the whole expression.
