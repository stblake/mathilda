### Worked examples

```mathematica
In[1]:= Clip[3.7]
Out[1]= 1
```

```mathematica
In[1]:= Clip[Pi]
Out[1]= 1
```

```mathematica
In[1]:= Clip[{-3, 0.5, 4}, {-1, 1}]
Out[1]= {-1, 0.5, 1}
```

```mathematica
In[1]:= Clip[15, {0, 10}, {-1, 1}]
Out[1]= 1
```

```mathematica
In[1]:= Clip[Infinity, {-2, 2}]
Out[1]= 2
```

### Notes

`Clip[x]` saturates `x` to the interval `[-1, 1]`; `Clip[x, {min, max}]`
saturates to `[min, max]`; and `Clip[x, {min, max}, {vmin, vmax}]` returns
`vmin`/`vmax` for out-of-range inputs, giving a piecewise ramp-and-saturate
profile. The first argument threads over lists, so a single call clips a whole
vector. Symbolic constants such as `Pi` are numericalized only to decide which
side of the interval they fall on, and `Infinity`/`-Infinity` clip to the upper
and lower replacement values respectively.
