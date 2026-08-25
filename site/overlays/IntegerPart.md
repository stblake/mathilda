---
status: Stable
---
### Worked examples

```mathematica
In[1]:= IntegerPart[2.7]
Out[1]= 2
```

```mathematica
In[1]:= IntegerPart[-2.7]
Out[1]= -2
```

```mathematica
In[1]:= IntegerPart[7/2]
Out[1]= 3
```

```mathematica
In[1]:= IntegerPart[{2.7, -2.7, 7/2, -7/2}]
Out[1]= {2, -2, 3, -3}
```

### Notes

`IntegerPart[x]` truncates **toward zero**, so `IntegerPart[-2.7] = -2` — unlike
`Floor`, which rounds toward `-Infinity` and would give `-3`. Exact inputs return
exact integers, and `IntegerPart` is `Listable`, threading over a vector of mixed
reals and rationals. Together with `FractionalPart` it splits any number as
`IntegerPart[x] + FractionalPart[x] == x`.
