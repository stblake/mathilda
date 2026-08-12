# SinIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SinIntegral[z]`**

gives the sine integral Si(z) = Integral\_0^z Sin\[t\]/t dt.

**`SinIntegral[+-Infinity] = +-Pi/2, SinIntegral[+-I Infinity] = +-I Infinity.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. SinIntegral\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[SinIntegral\[z\], z\] = Sinc\[z\]. Listable.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= SinIntegral[2.8]
Out[1]= 1.8321
```

```mathematica
In[1]:= N[SinIntegral[2], 50]
Out[1]= 1.6054129768026948485767201481985889408485834223285
```

```mathematica
In[1]:= SinIntegral[{-Infinity, Infinity, -I Infinity, I Infinity}]
Out[1]= {-1/2 Pi, 1/2 Pi, -I Infinity, I Infinity}
```

```mathematica
In[1]:= SinIntegral[2.5 + I]
Out[1]= 1.99549 + 0.222995 I
```

```mathematica
In[1]:= D[SinIntegral[x], x]
Out[1]= Sinc[x]
```

```mathematica
In[1]:= Series[SinIntegral[x], {x, 0, 7}]
Out[1]= x - 1/18 x^3 + 1/600 x^5 - 1/35280 x^7 + O[x]^8
```

```mathematica
In[1]:= Normal[Series[SinIntegral[x], {x, Infinity, 3}]]
Out[1]= 1/2 Pi - Sin[x]/x^2 + Cos[x] (-1/x + 2/x^3)
```

## Algorithm

```text
 Mathilda -- the sine integral  Si(z) = Int_0^z Sin[t]/t dt.

  SinIntegral[z]   Si(z)
```

Si is entire and odd, with no branch cuts. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-Pi/2, +-I Infinity, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (DLMF 6.6.5) is

```text
  Si(z) = Sum_{k>=0} (-1)^k z^(2k+1) / ((2k+1) (2k+1)!),
```

valid for all z. Its partial sums can reach magnitude ~e^|z| before the O(1)-sized answer emerges, so the MPFR paths add ~|z|/ln2 guard bits to absorb that cancellation exactly. For large |z| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 6.12.3)

```text
  Si(z) = Pi/2 - cos(z) f(z) - sin(z) g(z),
    f(z) ~ Sum (-1)^k (2k)!   / z^(2k+1) = 1/z - 2!/z^3 + ...,
    g(z) ~ Sum (-1)^k (2k+1)! / z^(2k+2) = 1/z^2 - 3!/z^4 + ...,
```

summed to the smallest term (optimal truncation). Both odd-symmetry reductions Si(-z) = -Si(z) fold negative real / left-half-plane inputs onto the right half-plane, keeping the asymptotic within its valid sector.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_nint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nint.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sinintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sinintegral.c)

## Notes & additional examples

### Notes

`SinIntegral[z]` is the sine integral `Si(z) = Integral_0^z Sin[t]/t dt`, an entire,
odd function with no branch cuts. Its derivative is [`Sinc`](Sinc.md), the cardinal
sine `Sin[z]/z`. On the imaginary axis `Si(I y) = I Shi(y)` in terms of the
hyperbolic sine integral, and as `x -> ±Infinity`, `Si(x) -> ±Pi/2`. A leading
negative is pulled out by odd symmetry (`SinIntegral[-x] = -SinIntegral[x]`). Numeric
evaluation uses a convergent Maclaurin series near the origin and an asymptotic
expansion for large `|z|`, at machine or arbitrary (MPFR) precision. Listable. See
also [`CosIntegral`](CosIntegral.md).
