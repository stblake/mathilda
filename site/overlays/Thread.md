### Worked examples

```mathematica
In[1]:= Thread[f[{a, b, c}]]
Out[1]= {f[a], f[b], f[c]}
```

```mathematica
In[1]:= Thread[{x, y, z} -> {1, 2, 3}]
Out[1]= {x -> 1, y -> 2, z -> 3}
```

```mathematica
In[1]:= Thread[f[{a, b}, {c, d}, x]]
Out[1]= {f[a, c, x], f[b, d, x]}
```

```mathematica
In[1]:= Thread[Equal[{a, b, c}, {1, 2, 3}]]
Out[1]= {a == 1, b == 2, c == 3}
```

### Notes

`Thread[f[args]]` distributes `f` over any lists in `args`, pairing them
positionally; non-list arguments are broadcast to every element. This is the
idiomatic way to build a rule list from parallel name/value lists
(`Thread[vars -> vals]`) or a system of equations from two vectors
(`Thread[Equal[lhs, rhs]]`). All list arguments must have the same length.
`Thread[f[args], h]` threads over a custom head `h` instead of `List`.
