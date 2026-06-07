### Worked examples

```mathematica
In[1]:= (Plus[##] &)[1, 2, 3]
Out[1]= 6

In[2]:= (f[##2] &)[a, b, c]
Out[2]= f[b, c]

In[3]:= FullForm[##]
Out[3]= SlotSequence[1]
```

### Notes

`##` (`SlotSequence[1]`) splices *all* arguments of the enclosing pure function into the surrounding expression; `##n` (Out[2]) splices the arguments from the n-th onward.
