### Worked examples

```mathematica
In[1]:= Inner[Times, {a, b}, {c, d}, Plus]
Out[1]= a c + b d
```

```mathematica
In[1]:= Inner[Times, {2, 3, 5}, {7, 11, 13}, Plus]
Out[1]= 112
```

```mathematica
In[1]:= Inner[f, {a, b}, {c, d}, g]
Out[1]= g[f[a, c], f[b, d]]
```

```mathematica
In[1]:= Inner[Times, {{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}, Plus]
Out[1]= {{19, 22}, {43, 50}}
```

### Notes

`Inner[f, list1, list2, g]` is the generalised dot product: `f` plays the role
of elementwise multiplication and `g` the role of summation. With `f = Times`
and `g = Plus` it reduces to ordinary `Dot`, so the matrix example above is just
the matrix product. Supplying symbolic `f` and `g` exposes the contraction
structure literally — `g[f[a, c], f[b, d]]` — which is useful for building
custom tensor operations (max-plus algebra, fuzzy logic, polynomial
convolutions, etc.). `Inner` contracts the last index of the first tensor with
the first index of the second.
