### Worked examples

```mathematica
In[1]:= FreeQ[x^2 + y^2, z]
Out[1]= True
```

```mathematica
In[1]:= FreeQ[x^2 + y^2, y]
Out[1]= False
```

```mathematica
In[1]:= FreeQ[D[Sin[x] Exp[x], x], Cos]
Out[1]= False
```

```mathematica
In[1]:= Select[Range[20], FreeQ[FactorInteger[#], {2, _}] &]
Out[1]= {1, 3, 5, 7, 9, 11, 13, 15, 17, 19}
```

### Notes

`FreeQ[expr, form]` tests whether *no* subexpression of `expr` matches the
pattern `form`, searching at every level. `form` may be a literal symbol or
a full pattern: the third example confirms that differentiating
`Sin[x] Exp[x]` introduces a `Cos` term. The last example is a structural
sieve — keeping only integers whose `FactorInteger` output is free of any
`{2, _}` pair recovers the odd numbers, with `FreeQ` matching against the
`{prime, exponent}` factor structure rather than a numeric value.
