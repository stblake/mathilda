---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= N[Sqrt[2]]
Out[1]= 1.41421
```

```mathematica
In[1]:= N[2/7, 15]
Out[1]= 0.2857142857142856
```

```mathematica
In[1]:= N[Pi, 40]
Out[1]= 3.1415926535897932384626433832795028841971
```

```mathematica
In[1]:= N[Zeta[3], 40]
Out[1]= 1.2020569031595942853997381615114499907651
```

```mathematica
In[1]:= N[Gamma[1/3], 35]
Out[1]= 2.67893853470774763365569294097467766
```

```mathematica
In[1]:= N[EulerGamma, 30]
Out[1]= 0.5772156649015328606065120900823
```

### Notes

`N[expr]` gives a machine-precision floating-point value, displayed to about six
significant digits. `N[expr, d]` requests approximately `d` digits of precision,
computed via arbitrary-precision arithmetic (so `N[Pi, 20]` returns the constant
to 20 digits). Exact inputs such as `Sqrt[2]`, `Pi`, `E`, and rationals are
converted to their numeric approximations. Note that machine-precision results
print at the default short width even when more digits are internally available.
