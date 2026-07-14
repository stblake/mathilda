### Worked examples

```mathematica
In[1]:= SinIntegral[2.8]
Out[1]= 1.8321
```

```mathematica
In[1]:= N[SinIntegral[2], 50]
Out[1]= 1.6054129768026948485767201481985889408485834223285
```

```mathematica
In[1]:= SinIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[1]= {-1/2 Pi, 1/2 Pi, -I Infinity, I Infinity}
```

```mathematica
In[1]:= SinIntegral[2.5 + I]
Out[1]= 1.99549 + 0.222995 I
```

```mathematica
In[1]:= D[SinIntegral[x], x]
Out[1]= Sinc[x]
```

```mathematica
In[1]:= Series[SinIntegral[x], {x, 0, 7}]
Out[1]= x - 1/18 x^3 + 1/600 x^5 - 1/35280 x^7 + O[x]^8
```

```mathematica
In[1]:= Normal[Series[SinIntegral[x], {x, Infinity, 3}]]
Out[1]= 1/2 Pi - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3)
```

### Notes

`SinIntegral[z]` is the sine integral `Si(z) = Integral_0^z Sin[t]/t dt`, an entire,
odd function with no branch cuts. Its derivative is [`Sinc`](Sinc.md), the cardinal
sine `Sin[z]/z`. On the imaginary axis `Si(I y) = I Shi(y)` in terms of the
hyperbolic sine integral, and as `x -> ±Infinity`, `Si(x) -> ±Pi/2`. A leading
negative is pulled out by odd symmetry (`SinIntegral[-x] = -SinIntegral[x]`). Numeric
evaluation uses a convergent Maclaurin series near the origin and an asymptotic
expansion for large `|z|`, at machine or arbitrary (MPFR) precision. Listable. See
also [`CosIntegral`](CosIntegral.md).
