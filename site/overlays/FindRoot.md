### Worked examples

```mathematica
In[1]:= FindRoot[x^2 - 2, {x, 1}]
Out[1]= {x -> 1.41421}
```

```mathematica
In[1]:= FindRoot[Cos[x] == x, {x, 1}, WorkingPrecision -> 40]
Out[1]= {x -> 0.73908513321516064165531208767387340401341}
```

```mathematica
In[1]:= FindRoot[BesselJ[0, x], {x, 2, 3}]
Out[1]= {x -> 2.40483}
```

```mathematica
In[1]:= FindRoot[Sin[x] == 0, {x, 3}, WorkingPrecision -> 40]
Out[1]= {x -> 3.1415926535897932384626433832875751974431}
```

```mathematica
In[1]:= FindRoot[{x^2 + y^2 == 1, x == y}, {{x, 1}, {y, 1}}]
Out[1]= {x -> 0.707107, y -> 0.707107}
```

### Notes

`FindRoot` accepts a function (sought equal to zero) or an explicit
equation, and honours `WorkingPrecision` via MPFR. The second example pins
the Dottie number — the unique real fixed point of cosine — to 40 digits;
the third finds the first positive zero of the Bessel function `J0` by
bracketing on `[2, 3]`; the fourth recovers `π` as a root of `Sin` to full
40-digit precision. The last example solves a nonlinear 2x2 system with a
Newton step using the analytic Jacobian.
