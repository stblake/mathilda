### Worked examples

```mathematica
In[1]:= StandardDeviation[{1, 2, 3, 4, 5}]
Out[1]= Sqrt[5/2]
```

Results are kept exact; here the (unbiased, n − 1) estimate of a small sample:

```mathematica
In[1]:= StandardDeviation[{2, 4, 4, 4, 5, 5, 7, 9}]
Out[1]= 4 Sqrt[2/7]
```

Numericalize the same estimate to 40 digits:

```mathematica
In[1]:= N[StandardDeviation[{2, 4, 4, 4, 5, 5, 7, 9}], 40]
Out[1]= 2.1380899352993950774764278470380281724321
```

It relates to `Variance` as its square root:

```mathematica
In[1]:= Variance[{1, 2, 3, 4, 5}]
Out[1]= 5/2
```

A constant sample has zero spread:

```mathematica
In[1]:= StandardDeviation[{1, 1, 1, 1}]
Out[1]= 0
```

### Notes

`StandardDeviation[data]` returns the sample (unbiased, divide-by-`n - 1`)
standard deviation, i.e. `Sqrt[Variance[data]]`. Exact inputs give exact
radical output, which `N[..., d]` evaluates to arbitrary precision. A list of
length 1 or a constant list yields `0`.
