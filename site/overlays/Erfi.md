### Worked examples

```mathematica
In[1]:= Erfi[0]
Out[1]= 0
```

```mathematica
In[1]:= Erfi[-x]
Out[1]= -Erfi[x]
```

```mathematica
In[1]:= N[Erfi[1], 30]
Out[1]= 1.650425758797542876025337729561
```

```mathematica
In[1]:= Series[Erfi[x], {x, 0, 7}]
Out[1]= 2/Sqrt[Pi] x + 2/3/Sqrt[Pi] x^3 + 1/5/Sqrt[Pi] x^5 + 1/21/Sqrt[Pi] x^7 + O[x]^8

In[2]:= D[Erfi[x], x]
Out[2]= (2 E^x^2)/Sqrt[Pi]
```

### Notes

`Erfi[z] = -I Erf[I z] = (2/Sqrt[Pi]) Integral_0^z e^(t^2) dt` is the imaginary
error function, an entire odd function with `Erfi[0] = 0`,
`Erfi[Infinity] = Infinity`, and `Erfi[I Infinity] = I`. Compared with `Erf`, the
sign of every Maclaurin coefficient is positive, reflecting the `+t^2` in the
integrand. Real and complex arguments evaluate numerically at machine or
arbitrary (MPFR) precision, the derivative is `(2/Sqrt[Pi]) E^(z^2)`, and `Erfi`
is `Listable`.
