### Worked examples

```mathematica
In[1]:= Erfc[0]
Out[1]= 1
```

```mathematica
In[1]:= N[Erfc[2], 40]
Out[1]= 0.0046777349810472658379307436327470713891081
```

```mathematica
In[1]:= N[Erfc[1 + I], 25]
Out[1]= -0.31615128169794764488027107 - 0.19045346923783468628410886*I
```

```mathematica
In[1]:= Series[Erfc[x], {x, 0, 5}]
Out[1]= 1 + -2/Sqrt[Pi] x + 2/3/Sqrt[Pi] x^3 + -1/5/Sqrt[Pi] x^5 + O[x]^6

In[2]:= D[Erfc[Sqrt[x]], x]
Out[2]= -E^(-x)/(Sqrt[Pi] Sqrt[x])
```

### Notes

`Erfc[z] = 1 - Erf[z]` is the complementary error function, with `Erfc[0] = 1`,
`Erfc[Infinity] = 0`, and `Erfc[-Infinity] = 2`. Real inputs evaluate through the
cancellation-free `erfc` of libm / MPFR — important in the right tail, where the
`N[Erfc[2], 40]` example keeps full precision instead of losing it to a `1 - Erf`
subtraction. Complex inputs route through `1 - Erf[z]` at machine or arbitrary
precision. The derivative is `D[Erfc[z], z] = -(2/Sqrt[Pi]) E^(-z^2)`, and `Erfc`
is `Listable`.
