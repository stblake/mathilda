### Worked examples

```mathematica
In[1]:= Decompose[x^4 + 6 x^2 + 10, x]
Out[1]= {10 + 6 x + x^2, x^2}

In[2]:= Decompose[x^6 + 6 x^4 + 12 x^2 + 8, x]
Out[2]= {8 + 12 x + 6 x^2 + x^3, x^2}
```

### Notes

`Decompose[poly, x]` returns `{p1, p2, ..., pk}` such that `poly == p1[p2[...[pk[x]]...]]`, with each `pi` of degree >= 2. Out[1] reads as `p1 = 10 + 6 x + x^2` composed with `p2 = x^2`. A polynomial with no nontrivial decomposition returns `{poly}` unchanged.
