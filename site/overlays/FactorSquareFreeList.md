### Worked examples

```mathematica
In[1]:= FactorSquareFreeList[(x - 1)^2 (x + 1)]  (* square-free factors with their multiplicities *)
```

```mathematica
In[2]:= FactorSquareFreeList[x^4 - 2 x^2 + 1]  (* (x^2 - 1)^2, grouped by multiplicity *)
```

```mathematica
In[3]:= FactorSquareFreeList[2 x^3 + 2 x^2]  (* the numerical factor leads *)
```

### Notes

`FactorSquareFreeList[poly]` gives the square-free factors of `poly` as `{factor,
exponent}` pairs, where the exponent is the multiplicity grouping repeated
factors — cheaper than full `FactorList` and enough when only multiplicities are
needed. It wraps `FactorSquareFree` (the Yun/Musser decomposition via GCDs of the
polynomial with its derivative) and forwards the `Extension` option verbatim. The
first pair is always the overall numerical factor `{c, 1}` (`{1, 1}` when there
is none). The square-free factors need not themselves be irreducible — here
`x^4 - 2 x^2 + 1` factors as `(x^2 - 1)^2`, not into `(x - 1)^2 (x + 1)^2`.
