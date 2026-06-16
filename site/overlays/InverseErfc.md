### Worked examples

```mathematica
In[1]:= InverseErfc[1]
Out[1]= 0
```

The endpoints of the domain `[0, 2]` are the infinities, with `InverseErfc[2] = -Infinity`:

```mathematica
In[1]:= InverseErfc[2]
Out[1]= -Infinity
```

High-precision evaluation, useful for tail quantiles where `Erf` underflows:

```mathematica
In[1]:= N[InverseErfc[1/1000000], 40]
Out[1]= 3.4589107372795000221509276359575695199155
```

The derivative is closed-form, `D[InverseErfc[z], z] == -(Sqrt[Pi]/2) E^(InverseErfc[z]^2)`:

```mathematica
In[1]:= D[InverseErfc[z], z]
Out[1]= -1/2 Sqrt[Pi] E^InverseErfc[z]^2
```

The reflection identity `InverseErf[1 - s] == InverseErfc[s]` is recognised symbolically:

```mathematica
In[1]:= InverseErf[1 - 3/10] == InverseErfc[3/10]
Out[1]= True
```

### Notes

`InverseErfc[s]` returns the `z` solving `Erfc[z] == s`, with
`InverseErfc[0] = Infinity`, `InverseErfc[1] = 0`, `InverseErfc[2] = -Infinity`.
Numerical values are produced only for real `s` in `[0, 2]`, at machine or
arbitrary (MPFR) precision. It is the natural function for accurate evaluation of
extreme normal-distribution quantiles, where `1 - Erf` would lose all precision
to cancellation.
