### Worked examples

The digamma function at positive integers reduces to an exact rational minus
Euler's constant:

```mathematica
In[1]:= PolyGamma[1]
Out[1]= -EulerGamma

In[2]:= PolyGamma[5]
Out[2]= 25/12 - EulerGamma
```

Higher derivatives `psi^(n)` at integer points give exact closed forms. The
trigamma at `1` is the Basel constant, and the odd-order values are rational
multiples of even powers of `Pi`:

```mathematica
In[1]:= PolyGamma[1, 1]
Out[1]= 1/6 Pi^2

In[2]:= PolyGamma[3, 1]
Out[2]= 1/15 Pi^4
```

The numeric paths cover arbitrary precision and complex arguments — here the
tetragamma value `psi^(2)(1+i)` to 30 digits:

```mathematica
In[1]:= N[PolyGamma[0, 3/2], 40]
Out[1]= 0.036489973978576520559023667001244432806843

In[2]:= N[PolyGamma[2, 1 + I], 30]
Out[2]= 0.3685529315879351717366345429807 + 0.7666528503450662124026953776316*I
```

The order `-1` is the special "integral of psi" case, the log-gamma function:

```mathematica
In[1]:= PolyGamma[-1, z]
Out[1]= LogGamma[z]
```

### Notes

`PolyGamma[z]` is the digamma function ψ(z), stored internally as
`PolyGamma[0, z]`; `PolyGamma[n, z]` is its n-th derivative ψ⁽ⁿ⁾(z). Positive
integer arguments reduce exactly: ψ(m) to a rational minus `EulerGamma`, and
ψ⁽ⁿ⁾(m) for odd `n` to a rational plus a rational multiple of `Pi^(n+1)` (even
orders are left symbolic). Non-positive integer arguments give `ComplexInfinity`
(the poles). Inexact real and complex arguments evaluate numerically at machine
or MPFR precision, and `PolyGamma[-1, z]` returns `LogGamma[z]`. Listable.
