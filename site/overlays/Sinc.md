### Worked examples

```mathematica
In[1]:= Sinc[0]
Out[1]= 1
```

```mathematica
In[1]:= Sinc[2.]
Out[1]= 0.454649
```

```mathematica
In[1]:= N[Sinc[2], 45]
Out[1]= 0.454648713412840847698009932955872421351127485
```

```mathematica
In[1]:= Sinc[1. + I]
Out[1]= 0.966711 - 0.331747 I
```

```mathematica
In[1]:= D[Sinc[x], x]
Out[1]= Cos[x]/x - Sin[x]/x^2
```

```mathematica
In[1]:= Series[Sinc[x], {x, 0, 6}]
Out[1]= 1 - x^2/6 + x^4/120 - x^6/5040 + O[x]^7
```

### Notes

`Sinc[z]` is the cardinal sine `Sin[z]/z`, with the removable singularity at the
origin filled in as `Sinc[0] = 1`. It is entire and even, and `Sinc[±Infinity] = 0`.
It appears as the derivative of the sine integral: `D[SinIntegral[z], z] = Sinc[z]`.
Numeric evaluation is at machine or arbitrary (MPFR) precision for both real and
complex arguments. Listable. See also [`SinIntegral`](SinIntegral.md).
