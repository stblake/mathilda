### Worked examples

```mathematica
In[1]:= UnitStep[-2]
Out[1]= 0

In[2]:= UnitStep[3]
Out[2]= 1
```

UnitStep is `Listable`, so it vectorizes over a table to produce a discrete step profile (here switching on at `k = 3`):

```mathematica
In[1]:= Table[UnitStep[k - 3], {k, 0, 6}]
Out[1]= {0, 0, 0, 1, 1, 1, 1}
```

Exact symbolic-real arguments are resolved by numerical certification — even transcendental comparisons collapse to an exact `0` or `1`:

```mathematica
In[1]:= UnitStep[Pi - 3]
Out[1]= 1

In[2]:= UnitStep[Log[2] - Log[3]]
Out[2]= 0
```

The multivariate form is the indicator of the nonnegative orthant, returning `1` only when no argument is negative:

```mathematica
In[1]:= UnitStep[1, -1, 2]
Out[1]= 0
```

### Notes

`UnitStep[x]` is `0` for `x < 0` and `1` for `x >= 0` (the value at `0` is `1`). The result is always exact: certifiable real arguments resolve numerically, while non-real or unresolved arguments are left unevaluated. UnitStep is `Listable` and `Orderless`.
