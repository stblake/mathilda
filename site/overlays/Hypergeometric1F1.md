### Worked examples

```mathematica
In[1]:= Hypergeometric1F1[1, 2, z]
Out[1]= (-1 + E^z)/z
```

A non-positive integer upper parameter terminates Kummer's series to a
polynomial (here a scaled Laguerre polynomial):

```mathematica
In[1]:= Hypergeometric1F1[-2, 1, z]
Out[1]= 1 - 2 z + 1/2 z^2
```

The numeric path agrees with the closed form to full precision — evaluating
`1F1[1, 2, 3]` and `(E^3 - 1)/3` to 40 digits:

```mathematica
In[1]:= N[Hypergeometric1F1[1, 2, 3], 40]
Out[1]= 6.3618456410625559136428432181939059656621

In[2]:= N[(E^3 - 1)/3, 40]
Out[2]= 6.3618456410625559136428432181939059656621
```

### Notes

`Hypergeometric1F1[a, b, z]` is Kummer's confluent hypergeometric function, equal to `HypergeometricPFQ[{a}, {b}, z]`, and converges for all `z`. A non-positive integer `a` truncates the series to a polynomial (the Laguerre/Hermite family); otherwise the function evaluates numerically at machine, MPFR, and complex precision.
