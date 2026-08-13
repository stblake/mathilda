# FresnelC

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FresnelC[z]`**

gives the Fresnel integral C(z) = Integral\_0^z Cos\[Pi t^2/2\] dt.

**`FresnelC[+-Infinity] = +-1/2, FresnelC[+-I Infinity] = +-I/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. FresnelC\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[FresnelC\[z\], z\] = Cos\[Pi z^2/2\]. Listable.

</details>

## Examples

_No verified examples yet for this function._

## Algorithm

Mathilda -- the Fresnel integrals (Pi/2-normalized).

```text
  FresnelC[z] = Int_0^z Cos[Pi t^2 / 2] dt
  FresnelS[z] = Int_0^z Sin[Pi t^2 / 2] dt
```

Both are entire and odd, with no branch cuts. FresnelC and FresnelS share one numeric kernel: the pair (C, S) is computed together and each builtin returns its component. Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-1/2, +-I/2, Indeterminate
  machine / arbitrary real ->  MPFR series (small |x|) or asymptotic (large)
  complex (any precision)  ->  the paired A/B ncpx series with guard bits
  everything else          ->  stays symbolic (return NULL)
```

Convergent Maclaurin series (valid for all z):

```text
  C(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m)   z^(4m+1) / ((2m)!   (4m+1)),
  S(z) = Sum_{m>=0} (-1)^m (Pi/2)^(2m+1) z^(4m+3) / ((2m+1)! (4m+3)).
```

Equivalently A(z) = C(z) + i S(z) = Sum_{k>=0} (i Pi/2)^k z^(2k+1)/(k!(2k+1)) and B(z) = C(z) - i S(z) is the same with i -> -i; then C = (A+B)/2 and S = (A-B)/(2i). The real path sums C and S as two real series directly; the complex path sums A and B together in one ncpx loop.

The partial sums reach magnitude ~e^((Pi/2)|z|^2) before the O(1)-sized answer emerges, so the MPFR paths add ~(Pi/2)|z|^2/ln2 guard bits to absorb that cancellation exactly. For large real |x| the convergent series is infeasible; there we use the asymptotic expansion (DLMF 7.12)

```text
  C(x) = 1/2 + f(x) sin(Pi x^2/2) - g(x) cos(Pi x^2/2),
  S(x) = 1/2 - f(x) cos(Pi x^2/2) - g(x) sin(Pi x^2/2),
    f(x) ~ (1/(Pi x))   Sum_j (-1)^j (4j-1)!! / (Pi x^2)^(2j),
    g(x) ~ (1/(Pi^2 x^3)) Sum_j (-1)^j (4j+1)!! / (Pi x^2)^(2j),
```

summed to the smallest term (optimal truncation). The asymptotic constant 1/2 is the value only in a sector around the real axis (Stokes phenomenon), so it is used for real inputs only; complex inputs always use the convergent series (correct everywhere, merely costlier for large |z|).

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
