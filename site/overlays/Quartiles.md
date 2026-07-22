### Worked examples

```mathematica
In[1]:= Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= {5/2, 9/2, 13/2}

In[2]:= Quartiles[{6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49}]
Out[2]= {81/4, 40, 171/4}
```

The estimates stay perfectly exact on rational data, so the interquartile range
(`q3 - q1`) of an arithmetic progression comes out in closed form:

```mathematica
In[1]:= q = Quartiles[Range[1, 1000]]
Out[1]= {501/2, 1001/2, 1501/2}

In[2]:= q[[3]] - q[[1]]
Out[2]= 500
```

An alternative quantile convention can be selected with a parameter list:

```mathematica
In[1]:= Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}, {{1/2, 0}, {0, 1}}]
Out[1]= {5/2, 9/2, 13/2}
```

### Notes

`Quartiles[data]` gives the three quartile estimates `{q1, q2, q3}` (lower,
median, upper) of the data. The default uses the standard quantile
definition; an optional parameter list `{{a, b}, {c, d}}` selects an
alternative quantile convention.
