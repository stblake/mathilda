# InverseErf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
InverseErf[s]
    gives the inverse error function: the z solving s = Erf[z].
InverseErf[z0, s]
    gives the inverse of the generalized error function Erf[z0, z].
InverseErf[0] = 0, InverseErf[1] = Infinity, InverseErf[-1] = -Infinity.
Odd in s. Numerical values are given only for real s in [-1, 1], at
machine or arbitrary (MPFR) precision; D[InverseErf[z], z] =
(Sqrt[Pi]/2) E^(InverseErf[z]^2). Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Exact special values: `InverseErf[0] = 0`, `InverseErf[1] = Infinity`,

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= InverseErf[0]
Out[1]= 0
```

High-precision evaluation for real arguments in `[-1, 1]`:

```mathematica
In[1]:= N[InverseErf[1/2], 40]
Out[1]= 0.47693627620446987338141835364313055980899
```

The Maclaurin series in powers of `Sqrt[Pi]`:

```mathematica
In[1]:= Series[InverseErf[x], {x, 0, 7}]
Out[1]= 1/2 Sqrt[Pi] x + 1/24 Pi^(3/2) x^3 + 7/960 Pi^(5/2) x^5 + 127/80640 Pi^(7/2) x^7 + O[x]^8
```

The derivative is closed-form, `D[InverseErf[z], z] == (Sqrt[Pi]/2) E^(InverseErf[z]^2)`:

```mathematica
In[1]:= D[InverseErf[z], z]
Out[1]= 1/2 Sqrt[Pi] E^InverseErf[z]^2
```

A statistical application: the two-sided 95% normal quantile is `Sqrt[2] InverseErf[2 p - 1]` with `p = 0.95`:

```mathematica
In[1]:= N[Sqrt[2] InverseErf[2 (95/100) - 1], 30]
Out[1]= 1.644853626951472714863848907989
```

### Notes

`InverseErf[s]` returns the `z` solving `Erf[z] == s`. It is odd in `s`, with
`InverseErf[0] = 0`, `InverseErf[1] = Infinity`, `InverseErf[-1] = -Infinity`.
Numerical values are produced only for real `s` in `[-1, 1]`, at machine or
arbitrary (MPFR) precision; symbolic arguments are returned unevaluated.
