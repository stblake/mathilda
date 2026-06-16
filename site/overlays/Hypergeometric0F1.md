### Worked examples

```mathematica
In[1]:= Hypergeometric0F1[1/2, z]
Out[1]= Cosh[2 Sqrt[z]]
```

With a different lower parameter the closed form switches to a hyperbolic
sine — `0F1` is the engine behind the Bessel and circular/hyperbolic
functions:

```mathematica
In[1]:= Hypergeometric0F1[3/2, z]
Out[1]= (1/2 Sinh[2 Sqrt[z]])/Sqrt[z]
```

Negative arguments give the trigonometric counterpart; with `z = -Pi^2/16`
the value of `0F1[3/2, z]` is `2/Pi` exactly, confirmed to 40 digits:

```mathematica
In[1]:= N[Hypergeometric0F1[3/2, -(Pi^2/16)]*Pi/2, 40]
Out[1]= 0.99999999999999999999999999999999999999991 + 2.3165577250480442772379354523138034341839e-41*I
```

### Notes

`Hypergeometric0F1[b, z]` is the confluent limit `HypergeometricPFQ[{}, {b}, z]`. It converges for all `z` and underlies the Bessel functions: `0F1[1/2, z] = Cosh[2 Sqrt[z]]` and `0F1[3/2, z] = Sinh[2 Sqrt[z]]/(2 Sqrt[z])`. The tiny imaginary residue in the last result is numerical noise from the radical of a negative argument.
