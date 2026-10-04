### Worked examples

```mathematica
In[1]:= SinhIntegral[0]
```

```mathematica
In[1]:= SinhIntegral[-2 x]  (* odd in its argument *)
```

```mathematica
In[1]:= SinhIntegral[Infinity]
```

```mathematica
In[1]:= D[SinhIntegral[x], x]  (* the integrand of the defining integral *)
```

```mathematica
In[1]:= Series[SinhIntegral[x], {x, 0, 7}]
```

```mathematica
In[1]:= N[SinhIntegral[1], 20]
```

```mathematica
In[1]:= N[SinhIntegral[2.0 + 1.0 I]]  (* a complex argument *)
```

### Notes

`SinhIntegral[z]` is the hyperbolic sine integral `Shi(z) = Int_0^z Sinh[t]/t
dt`, an entire, odd function and the imaginary-axis sibling of `SinIntegral`
(`Shi(z) = -I SinIntegral[I z]`). Its derivative is `Sinh[z]/z`, the integrand
of its defining integral.

Numeric evaluation sums the Maclaurin series for moderate `|z|` and switches to
an asymptotic expansion for large `|z|`; because every real-axis term is
positive there is no cancellation, and a real result that overflows a machine
double (`Shi` grows like `e^|z|`) is returned at extended exponent range. The
special values `SinhIntegral[±Infinity] = ±Infinity` and
`SinhIntegral[±I Infinity] = ±I Pi/2` are returned directly.
