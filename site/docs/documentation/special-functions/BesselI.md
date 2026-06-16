# BesselI

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BesselI[n, z]
    gives the modified Bessel function of the first kind I_n(z), the
    solution of z^2 y'' + z y' - (z^2 + n^2) y = 0 regular at the origin.
I_0(0) = 1, I_n(0) = 0 for integer n != 0; I_n grows like e^z as
z -> Inf and is even in n (I_{-n} = I_n). Has a branch cut along the
negative real z axis for non-integer n. Real and complex order and
argument evaluate numerically at machine or arbitrary (MPFR) precision;
D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BesselI[0, 2.0]
Out[1]= 2.27959

In[2]:= D[BesselI[n, x], x]
Out[2]= 1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
