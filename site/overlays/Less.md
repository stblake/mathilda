### Worked examples

```mathematica
In[1]:= 2 < 3
Out[1]= True

In[2]:= 5 < 5
Out[2]= False

In[3]:= Less[1, 10]
Out[3]= True
```

### Notes

`<` is the operator form of `Less`. A chained form such as `a < b < c` parses to `Inequality` and is decided pairwise.
