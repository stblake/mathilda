### Worked examples

```mathematica
In[1]:= OwnValues[a]
Out[1]= {}

In[2]:= a = 7
Out[2]= 7

In[3]:= OwnValues[a]
Out[3]= {7 -> 7}
```

### Notes

`OwnValues[s]` returns the direct value rules created by `s = ...`. An undefined symbol has an empty list; assigning `a = 7` stores a single own-value rule.
