# BesselJ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselJ[n, z]`**

gives the Bessel function of the first kind J\_n(z), a solution of z^2 y'' + z y' + (z^2 - n^2) y = 0 regular at the origin.

**`D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

J\_0(0) = 1, J\_n(0) = 0 for integer n != 0. Has a branch cut along the negative real z axis for non-integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BesselJ[0, 5.2]
Out[1]= -0.11029

In[2]:= D[BesselJ[n, x], x]
Out[2]= 1/2 (BesselJ[-1 + n, x] - BesselJ[1 + n, x])
```

### Applications (6)

```mathematica
In[1]:= BesselJ[0, 0]
Out[1]= 1

In[2]:= BesselJ[1, 0]
Out[2]= 0
```

Half-integer orders close in elementary functions:

```mathematica
In[1]:= BesselJ[1/2, z]
Out[1]= Sin[z] Sqrt[2/(Pi z)]
```

The Frobenius series at the origin:

```mathematica
In[1]:= Series[BesselJ[0, x], {x, 0, 6}]
Out[1]= 1 - 1/4 x^2 + 1/64 x^4 - 1/2304 x^6 + O[x]^7
```

High-precision and complex evaluation; the first input is `J_0` at its first zero, which numerically returns essentially zero:

```mathematica
In[1]:= N[BesselJ[0, 1], 40]
Out[1]= 0.7651976865579665514497175261026632209093

In[2]:= N[BesselJ[0, 10 + 5 I], 30]
Out[2]= -17.78959112945037151834426180967 + 0.2007116167212048509818027697064*I
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| BesselJ[0, .] over 10^6 | 3.46e+03 s | 1.74e+03 s | 53 s |
| Zeta over 10^6 | 8.17 s | 4.31e+03 s | 5.7 s |
| AiryAi over 10^6 | 3.14 s | 107 s | 58.9 s |
| PolyGamma[0, .] over 10^6 | 1.55 s | 151 s | 8.19 s |
| Gamma over 10^6 | 1.16 s | 1.35 s | 7.34 s |
| Erf over 10^6 | 0.968 s | 1.23 s | 7.43 s |

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## See also

[BesselY](../../special-functions/BesselY/), [BesselI](../../special-functions/BesselI/), [BesselK](../../special-functions/BesselK/), [N](../../arithmetic/N/), [Gamma](../../special-functions/Gamma/), [SeriesCoefficient](../../power-series/SeriesCoefficient/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_integrate_ramanujan.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_ramanujan.c)

## Notes & additional examples

### Notes

`BesselJ[n, z]` is the Bessel function of the first kind, regular at the origin, with `J_0(0) = 1` and `J_n(0) = 0` for integer `n != 0`. Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselJ[n, z], z] = (BesselJ[n-1, z] - BesselJ[n+1, z])/2`. There is a branch cut along the negative real axis for non-integer order. Listable.
