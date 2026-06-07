### Worked examples

```mathematica
In[1]:= ReplaceList[{a, b, c}, {Longest[x__], y__} :> {{x}, {y}}]
Out[1]= {{{a, b}, {c}}, {{a}, {b, c}}}

In[2]:= MatchQ[{1, 1}, {Longest[1 ..]}]
Out[2]= True
```

### Notes

`Longest[p]` forces a sequence pattern to consume as many elements as possible. In `ReplaceList` (Out[1]) the longest partition is tried first; contrast with `Shortest`, which orders the matches the other way.
