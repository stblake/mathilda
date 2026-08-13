# BesselI

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BesselI[n, z]`**

gives the modified Bessel function of the first kind I\_n(z), the solution of z^2 y'' + z y' - (z^2 + n^2) y = 0 regular at the origin.

**`D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2. Listable.`**

<details>
<summary>Notes</summary>

I\_0(0) = 1, I\_n(0) = 0 for integer n != 0; I\_n grows like e^z as z -\> Inf and is even in n (I\_{-n} = I\_n). Has a branch cut along the negative real z axis for non-integer n. Real and complex order and argument evaluate numerically at machine or arbitrary (MPFR) precision;

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BesselI[0, 2.0]
Out[1]= 2.27959

In[2]:= D[BesselI[n, x], x]
Out[2]= 1/2 (BesselI[-1 + n, x] + BesselI[1 + n, x])
```

### Applications (5)

```mathematica
In[3]:= BesselI[0, 0]
Out[3]= 1

In[4]:= BesselI[1/2, z]
Out[4]= Sinh[z] Sqrt[2/(Pi z)]

In[5]:= Series[BesselI[0, x], {x, 0, 6}]
Out[5]= 1 + 1/4 x^2 + 1/64 x^4 + 1/2304 x^6 + O[x]^7

In[6]:= N[BesselI[0, 1], 40]
Out[6]= 1.2660658777520083355982446252147175376077

In[7]:= N[BesselI[0, 2] BesselK[1, 2] + BesselI[1, 2] BesselK[0, 2], 30]
Out[7]= 0.5
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`, `ReadProtected`.

## References

**See also:** [BesselJ](../../special-functions/BesselJ/), [BesselK](../../special-functions/BesselK/), [N](../../arithmetic/N/), [Normal](../../power-series/Normal/), [Plus](../../arithmetic/Plus/), [Series](../../power-series/Series/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_besseli.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besseli.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`BesselI[n, z]` is the modified Bessel function of the first kind, regular at the origin, with `I_0(0) = 1` and `I_n(0) = 0` for integer `n != 0`. It grows like `e^z` and is even in `n` (`I_{-n} = I_n`). Real and complex order and argument evaluate at machine or MPFR precision; `D[BesselI[n, z], z] = (BesselI[n-1, z] + BesselI[n+1, z])/2`. Listable.
