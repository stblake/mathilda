### Worked examples

```mathematica
In[1]:= FactorInteger[60]
Out[1]= {{2, 2}, {3, 1}, {5, 1}}
```

```mathematica
In[1]:= FactorInteger[2^67 - 1]
Out[1]= {{193707721, 1}, {761838257287, 1}}
```

```mathematica
In[1]:= FactorInteger[20!]
Out[1]= {{2, 18}, {3, 8}, {5, 4}, {7, 2}, {11, 1}, {13, 1}, {17, 1}, {19, 1}}
```

```mathematica
In[1]:= FactorInteger[1000000000000066600000000000001]
Out[1]= {{1000000000000066600000000000001, 1}}
```

### Notes

`FactorInteger[n]` returns `{prime, exponent}` pairs. `2^67 - 1` reproduces F. N. Cole's celebrated 1903 factorisation of the Mersenne number M₆₇ into 193707721 × 761838257287, demonstrating the vendored GMP-ECM backend on a hard semiprime. The exponent vector of `20!` is exactly Legendre's formula applied to each prime ≤ 20. The thirty-one-digit "Belphegor prime" `10^30 + 666·10^14 + 1` is returned as a single factor, i.e. it is certified prime.
