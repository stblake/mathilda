# BesselY

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BesselY[n, z]
    gives the Bessel function of the second kind Y_n(z), the solution of
    z^2 y'' + z y' + (z^2 - n^2) y = 0 singular at the origin.
Y_0(0) = -Infinity, Y_n(0) = ComplexInfinity for integer n != 0; Y_n has
a logarithmic branch point at 0 and a branch cut along the negative real
z axis, with Y_{-n} = (-1)^n Y_n for integer n. Real and complex order
and argument evaluate numerically at machine or arbitrary (MPFR)
precision; D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2.
Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BesselY[0, 2.5]
Out[1]= 0.49807

In[2]:= D[BesselY[n, x], x]
Out[2]= 1/2 (BesselY[-1 + n, x] - BesselY[1 + n, x])
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[BesselY[1, 3.0]]
Out[1]= 0.324674
```

Half-integer orders close in elementary functions (the second-kind partner of `BesselJ[1/2, z]`):

```mathematica
In[1]:= BesselY[1/2, z]
Out[1]= -Cos[z] Sqrt[2/(Pi z)]
```

High-precision evaluation at the origin-neighbourhood and beyond:

```mathematica
In[1]:= N[BesselY[0, 1], 40]
Out[1]= 0.088256964215676957982926766023515162827815
```

The Wronskian of the two first-order solutions confirms `J_1(z) Y_0(z) - J_0(z) Y_1(z) = 2/(Pi z)`, here matching `2/(5 Pi)` at `z = 5`:

```mathematica
In[1]:= N[BesselJ[1, 5] BesselY[0, 5] - BesselJ[0, 5] BesselY[1, 5], 30]
Out[1]= 0.127323954473516268615107010698

In[2]:= N[2/(5 Pi), 30]
Out[2]= 0.127323954473516268615107010698
```

### Notes

`BesselY[n, z]` is the Bessel function of the second kind, singular at the origin: `Y_0(0) = -Infinity`, `Y_n(0) = ComplexInfinity` for integer `n != 0`, with a logarithmic branch point at 0 and a branch cut along the negative real axis (`Y_{-n} = (-1)^n Y_n` for integer `n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselY[n, z], z] = (BesselY[n-1, z] - BesselY[n+1, z])/2`. Listable.
