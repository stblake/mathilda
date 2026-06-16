---
status: Stable
references:
  - "Harold Abelson and Gerald Jay Sussman, *Structure and Interpretation of Computer Programs*, 2nd ed., §2.2.1 (sequence mapping)."
---
### Worked examples

```mathematica
In[1]:= Map[f, {a, b, c}]
Out[1]= {f[a], f[b], f[c]}
```

```mathematica
In[1]:= #^2 & /@ {1, 2, 3, 4}
Out[1]= {1, 4, 9, 16}
```

```mathematica
In[1]:= Map[Reverse, {{1, 2}, {3, 4}}]
Out[1]= {{2, 1}, {4, 3}}
```

```mathematica
In[1]:= Map[f, {{a}, {b}}, {2}]
Out[1]= {{f[a]}, {f[b]}}
```

```mathematica
In[1]:= Map[Total, {{1, 2, 3}, {4, 5, 6}}]
Out[1]= {6, 15}
```

```mathematica
In[1]:= Map[#^2 &, x + y + z]
Out[1]= x^2 + y^2 + z^2
```

### Notes

`f /@ expr` is the operator shorthand for `Map[f, expr]`. By default the function
is applied at level 1, i.e. to the immediate elements; a level specification such
as `{2}` reaches deeper into nested lists. `Map` works on any expression, not
only `List` — the head is preserved while each argument is wrapped by `f`. Pure
functions (`#^2 &`) are the idiomatic first argument.
