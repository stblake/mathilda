### Worked examples

```mathematica
In[1]:= 1 < 2 < 3
Out[1]= True

In[2]:= 1 < 5 < 3
Out[2]= False

In[3]:= 2 <= 2 < 5
Out[3]= True
```

### Notes

A chained comparison such as `a < b <= c` parses to the canonical `Inequality` form and holds only when every adjacent pair holds. Mixed operators (e.g. `<=` and `<`) may be combined in a single chain.
