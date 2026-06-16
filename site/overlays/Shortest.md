### Worked examples

With `/.` only the single first match is taken; `Shortest` makes that the
fewest-element partition, while `Longest` makes it the most:

```mathematica
In[1]:= {a, b, c} /. {Shortest[x__], y__} :> {{x}, {y}}
Out[1]= {{a}, {b, c}}

In[2]:= {a, b, c} /. {Longest[x__], y__} :> {{x}, {y}}
Out[2]= {{a, b}, {c}}
```

`ReplaceList` enumerates every match, revealing the order in which partitions are
attempted — `Shortest` tries the smallest first, `Longest` the largest:

```mathematica
In[1]:= ReplaceList[{a, b, c}, {Shortest[x__], y__} :> {{x}, {y}}]
Out[1]= {{{a}, {b, c}}, {{a, b}, {c}}}

In[2]:= ReplaceList[{a, b, c}, {Longest[x__], y__} :> {{x}, {y}}]
Out[2]= {{{a, b}, {c}}, {{a}, {b, c}}}
```

### Notes

`Shortest[p]` forces a sequence pattern to consume as few elements as possible, so the shortest partition is tried first (Out[1]). It is the dual of `Longest` (Out[2]).
