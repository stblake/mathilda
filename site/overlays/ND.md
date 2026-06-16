### Worked examples

```mathematica
In[1]:= ND[Sin[x], x, 1]
Out[1]= 0.540302
```

```mathematica
In[1]:= ND[Gamma[x], x, 1]
Out[1]= -0.577216
```

```mathematica
In[1]:= ND[BesselJ[0, x], x, 2]
Out[1]= -0.576725
```

```mathematica
In[1]:= ND[Tan[x], {x, 2}, 1]
Out[1]= 11.4484
```

### Notes

`ND[expr, x, x0]` numerically differentiates `expr` at `x = x0`. The first case
recovers `Cos[1] = 0.540302`. The Gamma example gives `Gamma'[1] = -EulerGamma`,
since `PolyGamma[0, 1] = -EulerGamma`. The Bessel example uses the identity
`BesselJ[0, x]' = -BesselJ[1, x]`, so the value is `-BesselJ[1, 2]`. The
`{x, 2}` form takes the second derivative. The default `Method -> EulerSum`
applies Richardson extrapolation to finite differences; `Method -> NIntegrate`
uses Cauchy's integral formula and allows fractional or complex orders. `ND`
cannot recognise small numbers that should be zero — `Chop` if needed.
