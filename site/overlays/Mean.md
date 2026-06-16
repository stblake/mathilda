### Worked examples

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2
```

```mathematica
In[1]:= Mean[{a, b, c}]
Out[1]= 1/3 (a + b + c)
```

```mathematica
In[1]:= Mean[{1/2, 1/3, 1/6}]
Out[1]= 1/3
```

```mathematica
In[1]:= Mean[Table[k^2, {k, 1, 10}]]
Out[1]= 77/2
```

### Notes

`Mean[data]` is the arithmetic mean — the sum of the elements divided by their
count. It works symbolically as well as numerically: `Mean[{a, b, c}]` returns
the exact closed form `(a + b + c)/3`. Numeric data stays in exact rational
arithmetic, so `Mean[{1, 2, 3, 4}]` is `5/2` (not `2.5`) and the mean of the
first ten squares is `77/2`, with no round-off. Combined with generators like
`Table` and `Range`, `Mean` gives exact averages of structured data sets.
