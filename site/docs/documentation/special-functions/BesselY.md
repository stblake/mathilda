# BesselY

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselY[n, z]`**

gives the Bessel function of the second kind Y\_n(z), the solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.

<details>
<summary>Notes</summary>

Y\_0(0) = -Infinity, Y\_n(0) = ComplexInfinity for integer n != 0; Y\_n has a logarithmic branch point at 0 and a branch cut along the negative real z axis, with Y\_{-n} = (-1)^n Y\_n for integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision; D\[BesselY\[n, z\], z\] = (BesselY\[n-1, z\] - BesselY\[n+1, z\])/2. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])
```

### Applications (5)

```mathematica
In[3]:= N[BesselY[1, 3.0]]
Out[3]= 0.324674

In[4]:= BesselY[1/2, z]
Out[4]= -Cos[z] Sqrt[2/(Pi z)]

In[5]:= N[BesselY[0, 1], 40]
Out[5]= 0.088256964215676957982926766023515162827815

In[6]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5], 30]
Out[6]= 0.127323954473516268615107010698

In[7]:= N[2/(5 Pi), 30]
Out[7]= 0.127323954473516268615107010698
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [BesselK](../../special-functions/BesselK/), [BesselI](../../special-functions/BesselI/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_limit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_limit.c)

## Notes & additional examples

### Notes

`BesselY[n, z]` is the Bessel function of the second kind, singular at the origin: `Y_0(0) = -Infinity`, `Y_n(0) = ComplexInfinity` for integer `n != 0`, with a logarithmic branch point at 0 and a branch cut along the negative real axis (`Y_{-n} = (-1)^n Y_n` for integer `n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`. Listable.
