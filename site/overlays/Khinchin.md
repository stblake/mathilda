### Worked examples

```mathematica
In[1]:= N[Khinchin]
Out[1]= 2.68545
```

Khinchin's constant evaluates to arbitrary precision via its convergent product over partial quotients:

```mathematica
In[1]:= N[Khinchin, 60]
Out[1]= 2.685452001065306445309714835481795693820382293994462953051151
```

It is a true symbolic constant — `NumericQ` is `True` and its derivative vanishes:

```mathematica
In[1]:= NumericQ[Khinchin]
Out[1]= True

In[2]:= D[Khinchin, x]
Out[2]= 0
```

### Notes

`Khinchin` is Khinchin's (Khintchine's) constant `K ~= 2.68545`, the limiting
geometric mean of the partial quotients in the continued-fraction expansion of
almost every real number: `K = Product[(1 + 1/(s (s + 2)))^Log2[s], {s, 1,
Infinity}]`. It carries the `Constant` and `Protected` attributes, so it stays
symbolic until `N[Khinchin, prec]` evaluates it to the requested precision.
