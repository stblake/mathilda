### Worked examples

```mathematica
In[1]:= AiryBi[0]
Out[1]= 1/(3^(1/6) Gamma[2/3])
```

```mathematica
In[1]:= N[AiryBi[0], 40]
Out[1]= 0.6149266274460007351509223690936135535947
```

```mathematica
In[1]:= D[AiryBi[z], z]
Out[1]= AiryBiPrime[z]
```

```mathematica
In[1]:= N[AiryBi[2.0 + 1.0 I], 20]
Out[1]= 0.778230383757041677129 + 2.50509630006410244363*I
```

### Notes

`AiryBi[z]` is the dominant solution of the Airy equation `y'' == z y`, growing
exponentially as `z -> +Infinity` while `AiryBi[-Infinity] == 0`. Its exact
value at the origin is `1/(3^(1/6) Gamma[2/3])`, and `D[AiryBi[z], z]` returns
`AiryBiPrime[z]`. Complex arguments are evaluated to the requested MPFR
precision; `AiryBi` is Listable.
