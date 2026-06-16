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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[BesselK[1, 2.5]]
Out[1]= 0.0738908
```

Half-integer orders give exponentially decaying elementary forms:

```mathematica
In[1]:= BesselK[1/2, z]
Out[1]= E^(-z) Sqrt[(1/2 Pi)/z]
```

High-precision real and complex evaluation:

```mathematica
In[1]:= N[BesselK[0, 1], 40]
Out[1]= 0.42102443824070833333562737921260903613623

In[2]:= N[BesselK[0, 3 + I], 30]
Out[2]= 0.01383067506051671850189255536523 - 0.03098977854031822729465893945728*I
```

The Wronskian with `BesselI` confirms `I_0(z) K_1(z) + I_1(z) K_0(z) = 1/z`, here `1/2` at `z = 2`:

```mathematica
In[1]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]
Out[1]= 0.5
```

### Notes

`BesselK[n, z]` is the modified Bessel function of the second kind, decaying like `e^{-z}` and even in `n` (`K_{-n} = K_n`). `K_0(0) = Infinity`, `K_n(0) = ComplexInfinity`, with a branch cut along the negative real axis. Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselK[n, z], z] = -(BesselK[n-1, z] + BesselK[n+1, z])/2`. Listable.
