---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= MatchQ[3, _Integer]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[f[a, b], f[_, _]]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[x^2, _^_]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[{1, 2, 3}, {__Integer}]
Out[1]= True
```

### Notes

`MatchQ[expr, form]` tests whether `expr` matches the pattern `form`, returning
`True` or `False`. It supports the full pattern language: typed blanks
(`_Integer`), structural patterns (`f[_, _]`, `_^_`), and sequence variables
(`__Integer` for one-or-more integers, `___` for zero-or-more). Because `x^2` is
stored as `Power[x, 2]`, it matches the structural pattern `_^_`. `MatchQ` is the
predicate underlying filters like `Cases` and conditional rules.
