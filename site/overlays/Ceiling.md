### Worked examples

```mathematica
In[1]:= Ceiling[7/2]
Out[1]= 4

In[2]:= Ceiling[-2.7]
Out[2]= -2

In[3]:= Ceiling[17, 5]
Out[3]= 20
```

```mathematica
In[1]:= Ceiling[Pi]
Out[1]= 4
```

```mathematica
In[1]:= Ceiling[E^2]
Out[1]= 8
```

```mathematica
In[1]:= Ceiling[100/7, 5]
Out[1]= 15
```

### Notes

`Ceiling[x]` rounds toward `+Infinity`; the two-argument `Ceiling[x, a]` gives the smallest multiple of `a` at least `x`. Exact inputs return exact integers. Symbolic constants are resolved exactly: `Ceiling[Pi]` is `4` and `Ceiling[E^2]` is `8` (since `e^2 ≈ 7.389`). The two-argument form rounds up to a lattice point, so `Ceiling[100/7, 5]` (with `100/7 ≈ 14.29`) gives the next multiple of 5, namely `15`.
