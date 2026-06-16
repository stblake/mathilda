### Worked examples

```mathematica
In[1]:= Array[f, 5]
Out[1]= {f[1], f[2], f[3], f[4], f[5]}

In[2]:= Array[#^2 &, 4]
Out[2]= {1, 4, 9, 16}

In[3]:= Array[f, 3, 0]
Out[3]= {f[0], f[1], f[2]}
```

A two-argument function over a dimension spec builds a matrix; here the identity matrix as a Kronecker delta:

```mathematica
In[1]:= Array[Boole[#1 == #2] &, {3, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

Building the 4x4 Hilbert matrix and taking its determinant gives the famously tiny rational, evidence of its near-singularity:

```mathematica
In[1]:= Det[Array[1/(#1 + #2 - 1) &, {4, 4}]]
Out[1]= 1/6048000
```

### Notes

`Array[f, n]` builds a length-`n` list by applying `f` to indices `1..n`; an optional third argument sets the starting index. A dimension list `{n1, n2, ...}` produces a nested array with elements `f[i1, i2, ...]`.
