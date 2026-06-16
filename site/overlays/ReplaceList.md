### Worked examples

```mathematica
In[1]:= ReplaceList[{a, b, c, d}, {x___, y___} :> {{x}, {y}}]
Out[1]= {{{}, {a, b, c, d}}, {{a}, {b, c, d}}, {{a, b}, {c, d}}, {{a, b, c}, {d}}, {{a, b, c, d}, {}}}

In[2]:= ReplaceList[a + b + c, x_ + y_ :> {x, y}]
Out[2]= {{a, b + c}, {b, a + c}, {c, a + b}, {a + b, c}, {a + c, b}, {b + c, a}}

In[3]:= ReplaceList[{a, b, c, d}, {x___, y___} :> {{x}, {y}}, 2]
Out[3]= {{{}, {a, b, c, d}}, {{a}, {b, c, d}}}
```

```mathematica
In[1]:= ReplaceList[{1, 2, 3, 4}, {x___, y_, z_, w___} /; y + z == 5 :> {y, z}]
Out[1]= {{2, 3}}
```

```mathematica
In[1]:= ReplaceList[{1, 2, 3, 4, 5}, {a___, b_, c___} /; b == 3 :> {{a}, {c}}]
Out[1]= {{{1, 2}, {4, 5}}}
```

### Notes

Unlike `Replace`, which returns only the first match, `ReplaceList` enumerates *every* way the rule can match — useful for sequence patterns (`__`, `___`) that admit multiple partitions. A third argument caps the number of results returned.
