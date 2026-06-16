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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= BesselI[0, 0]
Out[1]= 1
```

Half-integer orders close in hyperbolic functions:

```mathematica
In[1]:= BesselI[1/2, z]
Out[1]= Sinh[z] Sqrt[2/(Pi z)]
```

The all-positive Frobenius series at the origin (contrast with the alternating series of `BesselJ`):

```mathematica
In[1]:= Series[BesselI[0, x], {x, 0, 6}]
Out[1]= 1 + 1/4 x^2 + 1/64 x^4 + 1/2304 x^6 + O[x]^7
```

High-precision evaluation, and a Wronskian identity with `BesselK`: `I_0(z) K_1(z) + I_1(z) K_0(z) = 1/z`, which at `z = 2` gives exactly `1/2`:

```mathematica
In[1]:= N[BesselI[0, 1], 40]
Out[1]= 1.2660658777520083355982446252147175376077

In[2]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]
Out[2]= 0.5
```

### Notes

`BesselI[n, z]` is the modified Bessel function of the first kind, regular at the origin, with `I_0(0) = 1` and `I_n(0) = 0` for integer `n != 0`. It grows like `e^z` and is even in `n` (`I_{-n} = I_n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2`. Listable.
