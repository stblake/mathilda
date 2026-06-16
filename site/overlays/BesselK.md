### Worked examples

```mathematica
In[1]:= N[BesselK[1, 2.5]]
Out[1]= 0.0738908
```

Half-integer orders give exponentially decaying elementary forms:

```mathematica
In[1]:= BesselK[1/2, z]
Out[1]= E^(-z) Sqrt[(1/2 Pi)/z]
```

High-precision real and complex evaluation:

```mathematica
In[1]:= N[BesselK[0, 1], 40]
Out[1]= 0.42102443824070833333562737921260903613623

In[2]:= N[BesselK[0, 3 + I], 30]
Out[2]= 0.01383067506051671850189255536523 - 0.03098977854031822729465893945728*I
```

The Wronskian with `BesselI` confirms `I_0(z) K_1(z) + I_1(z) K_0(z) = 1/z`, here `1/2` at `z = 2`:

```mathematica
In[1]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]
Out[1]= 0.5
```

### Notes

`BesselK[n, z]` is the modified Bessel function of the second kind, decaying like `e^{-z}` and even in `n` (`K_{-n} = K_n`). `K_0(0) = Infinity`, `K_n(0) = ComplexInfinity`, with a branch cut along the negative real axis. Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2`. Listable.
