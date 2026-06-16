### Worked examples

```mathematica
In[1]:= ExpToTrig[Exp[x]]
Out[1]= Cosh[x] + Sinh[x]
```

```mathematica
In[1]:= ExpToTrig[Exp[I x]]
Out[1]= Cos[x] + I Sin[x]
```

```mathematica
In[1]:= ExpToTrig[(E^x - E^(-x))/2]
Out[1]= Sinh[x]
```

```mathematica
In[1]:= ExpToTrig[E^(I x) + E^(-I x)]
Out[1]= 2 Cos[x]
```

### Notes

`ExpToTrig` rewrites exponentials in terms of circular and hyperbolic
functions: a real exponent yields `Cosh + Sinh`, an imaginary one yields
`Cos + I Sin`. It recognises the classical combinations, folding the
half-difference `(E^x - E^(-x))/2` back to `Sinh[x]` and the sum
`E^(I x) + E^(-I x)` to `2 Cos[x]`.
