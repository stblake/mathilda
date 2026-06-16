### Worked examples

```mathematica
In[1]:= N[BesselY[1, 3.0]]
Out[1]= 0.324674
```

Half-integer orders close in elementary functions (the second-kind partner of `BesselJ[1/2, z]`):

```mathematica
In[1]:= BesselY[1/2, z]
Out[1]= -Cos[z] Sqrt[2/(Pi z)]
```

High-precision evaluation at the origin-neighbourhood and beyond:

```mathematica
In[1]:= N[BesselY[0, 1], 40]
Out[1]= 0.088256964215676957982926766023515162827815
```

The Wronskian of the two first-order solutions confirms `J_1(z) Y_0(z) - J_0(z) Y_1(z) = 2/(Pi z)`, here matching `2/(5 Pi)` at `z = 5`:

```mathematica
In[1]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5], 30]
Out[1]= 0.127323954473516268615107010698

In[2]:= N[2/(5 Pi), 30]
Out[2]= 0.127323954473516268615107010698
```

### Notes

`BesselY[n, z]` is the Bessel function of the second kind, singular at the origin: `Y_0(0) = -Infinity`, `Y_n(0) = ComplexInfinity` for integer `n != 0`, with a logarithmic branch point at 0 and a branch cut along the negative real axis (`Y_{-n} = (-1)^n Y_n` for integer `n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`. Listable.
