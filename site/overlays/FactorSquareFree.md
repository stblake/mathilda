### Worked examples

```mathematica
In[1]:= FactorSquareFree[x^5 - x^4 - x + 1]
Out[1]= (-1 + x)^2 (1 + x + x^2 + x^3)
```

```mathematica
In[1]:= FactorSquareFree[(x^2+1)^3 (x-1)^2]
Out[1]= (-1 + x)^2 (1 + x^2)^3
```

```mathematica
In[1]:= FactorSquareFree[x^8 + 4 x^6 + 6 x^4 + 4 x^2 + 1]
Out[1]= (1 + x^2)^4
```

### Notes

`FactorSquareFree` groups a polynomial into pairwise-coprime square-free factors carrying their multiplicities, using the Yun / Musser decomposition (GCDs of the polynomial with its derivative). It does not split factors that are square-free but reducible: `1 + x + x^2 + x^3` is left intact in the first example even though it factors as `(1+x)(1+x^2)`. The last example recognises `(1 + x^2)^4` as a perfect fourth power without ever calling full `Factor`, which is the point — it is cheaper than factorisation when only multiplicities are needed.
