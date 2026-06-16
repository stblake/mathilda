# BesselK

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BesselK[n, z]
    gives the modified Bessel function of the second kind K_n(z), a
    solution of z^2 y'' + z y' - (z^2 + n^2) y = 0.
K_n(z) is even in n (K_{-n} = K_n) and decays like e^{-z} as z -> Inf.
K_0(0) = Infinity, K_n(0) = ComplexInfinity. Has a branch cut along the
negative real z axis. Real and complex order and argument evaluate
numerically at machine or arbitrary (MPFR) precision;
D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2. Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= BesselK[0, 4.0]
Out[1]= 0.0111597

In[2]:= D[BesselK[n, x], x]
Out[2]= -1/2 (BesselK[-1 + n, x] + BesselK[1 + n, x])
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
