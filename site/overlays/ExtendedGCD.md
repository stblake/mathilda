### Worked examples

```mathematica
In[1]:= ExtendedGCD[12, 18]
Out[1]= {6, {-1, 1}}

In[2]:= ExtendedGCD[15, 25, 35]
Out[2]= {5, {2, -1, 0}}
```

### Notes

`ExtendedGCD[n1, ...]` returns `{g, {r1, r2, ...}}` where `g == GCD[n1, ...]` and `g == r1 n1 + r2 n2 + ...`. For Out[1], `6 == (-1)(12) + (1)(18)`; the multi-argument form folds `mpz_gcdext` pairwise.
