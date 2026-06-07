### Worked examples

```mathematica
In[1]:= Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= {5/2, 9/2, 13/2}

In[2]:= Quartiles[{6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49}]
Out[2]= {81/4, 40, 171/4}
```

### Notes

`Quartiles[data]` gives the three quartile estimates `{q1, q2, q3}` (lower,
median, upper) of the data. The default uses Mathematica's standard quantile
definition; an optional parameter list `{{a, b}, {c, d}}` selects an
alternative quantile convention.
