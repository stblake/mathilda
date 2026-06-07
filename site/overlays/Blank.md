### Worked examples

```mathematica
In[1]:= MatchQ[5, _]
Out[1]= True

In[2]:= MatchQ[g[a], _g]
Out[2]= True

In[3]:= FullForm[x_]
Out[3]= Pattern[x, Blank[]]
```

### Notes

`_` (`Blank[]`) matches any single expression; `_h` (`Blank[h]`) matches any single expression whose head is `h` (Out[2]). A named blank `x_` is `Pattern[x, Blank[]]` (Out[3]), binding the match to `x`.
