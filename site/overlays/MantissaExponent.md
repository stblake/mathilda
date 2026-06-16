### Worked examples

```mathematica
In[1]:= MantissaExponent[123.45]
Out[1]= {0.12345, 3}
```

```mathematica
In[1]:= MantissaExponent[7/3]
Out[1]= {7/30, 1}
```

```mathematica
In[1]:= MantissaExponent[1024, 2]
Out[1]= {1/2, 11}
```

```mathematica
In[1]:= MantissaExponent[N[Pi, 30]]
Out[1]= {0.3141592653589793238462643383278, 1}
```

```mathematica
In[1]:= MantissaExponent[255, 16]
Out[1]= {255/256, 2}
```

### Notes

`MantissaExponent[x]` returns `{m, e}` with `x = m * 10^e` and `1/10 <= |m| < 1` (or `{0, 0}` when `x` is 0). `MantissaExponent[x, b]` uses base `b`, so `1/b <= |m| < 1`. Exact inputs keep an exact `Rational` mantissa (`7/3 -> {7/30, 1}`); inexact inputs keep their full working precision (the `N[Pi, 30]` mantissa carries 30 digits). Only integer bases `>= 2` are currently supported.
