### Worked examples

```mathematica
In[1]:= DeleteCases[{1,2,3,2},2]
Out[1]= {1, 3}
```

```mathematica
In[1]:= DeleteCases[{1,a,2,b,3},_Integer]
Out[1]= {a, b}
```

```mathematica
In[1]:= DeleteCases[{{1,2},{a,b},{3,c}}, {_Integer, _Integer}]
Out[1]= {{a, b}, {3, c}}
```

```mathematica
In[1]:= DeleteCases[{a + b, c d, e^2, f}, _Power, Infinity]
Out[1]= {a + b, c d, f}
```

### Notes

Removes every element matching the pattern; combine with `_Integer` or `_?OddQ` to filter by head or predicate. The pattern can be structural — `{_Integer, _Integer}` deletes only the sublists that are pairs of integers. With a level specification such as `Infinity`, `DeleteCases` descends into subexpressions and removes matching parts at any depth (here every `Power` subterm), traversing in depth-first post-order. The default level is `{1}` (top-level elements).
