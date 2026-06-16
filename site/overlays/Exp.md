### Worked examples

```mathematica
In[1]:= Exp[0]
Out[1]= 1

In[2]:= Exp[I Pi]
Out[2]= -1

In[3]:= Exp[Log[x]]
Out[3]= x
```

```mathematica
In[1]:= Exp[2 Log[x]]
Out[1]= x^2
```

```mathematica
In[1]:= D[Exp[Sin[x]], x]
Out[1]= Cos[x] E^Sin[x]
```

```mathematica
In[1]:= Series[Exp[x], {x, 0, 6}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + 1/720 x^6 + O[x]^7
```

```mathematica
In[1]:= N[Exp[1], 50]
Out[1]= 2.71828182845904523536028747135266249775724709369996
```

### Notes

`Exp[z]` is `E^z`; it inverts `Log` and evaluates Euler's identity exactly.
Logarithmic arguments collapse (`Exp[2 Log[x]] -> x^2`), it differentiates and
series-expands symbolically, and numeric inputs route to libm / MPFR for
arbitrary precision. Exp is Listable.
