### Worked examples

```mathematica
In[1]:= Min[3, 7, 2]
Out[1]= 2
```

```mathematica
In[1]:= Min[1/3, 2/7, 5/11]
Out[1]= 2/7
```

```mathematica
In[1]:= Min[x, 0, Infinity]
Out[1]= Min[0, x]
```

```mathematica
In[1]:= Min[{}]
Out[1]= Infinity
```

### Notes

`Min[x1, x2, ...]` returns the numerically smallest argument, and `Min` of
several lists returns the smallest element across all of them. Comparisons are
exact, so rationals are ordered without converting to floating point —
`Min[1/3, 2/7, 5/11]` correctly picks `2/7`. With symbolic arguments `Min` stays
unevaluated but still prunes operands it can decide: `Min[x, 0, Infinity]` drops
`Infinity` (which can never be a minimum) and returns `Min[0, x]`. The empty case
`Min[{}]` returns `Infinity`, the identity element of minimisation — the value
that leaves any subsequent `Min` unchanged.
