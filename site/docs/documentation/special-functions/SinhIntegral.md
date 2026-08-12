# SinhIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SinhIntegral[z]`**

gives the hyperbolic sine integral Shi(z) = Integral\_0^z Sinh\[t\]/t dt.

**`SinhIntegral[+-Infinity] = +-Infinity, SinhIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

An entire, odd function with no branch cuts. SinhIntegral\[0\] = 0, Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[SinhIntegral\[z\], z\] = Sinh\[z\]/z. Listable.

</details>

## Examples

_No verified examples yet for this function._

## Algorithm

```text
 Mathilda -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.

  SinhIntegral[z]   Shi(z)
```

Shi is entire and odd, with no branch cuts. It is the imaginary-axis sibling of the sine integral: Shi(z) = -i Si(i z). Evaluation is layered so each kind of argument takes the cheapest route:

```text
  exact special values     ->  0, +-Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Shi(z) = Sum_{k>=0} z^(2k+1) / ((2k+1) (2k+1)!),
```

valid for all z. Unlike Si, every term is positive, so on the real axis there is NO catastrophic cancellation -- the partial sums climb monotonically to the O(e^|z|)-sized answer and only a small fixed guard is needed. (Complex inputs still cancel like Si, so the complex paths keep the ~|z|/ln2 guard.) For large

```text
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Shi(z) ~ cosh(z) F(z) + sinh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

summed to the smallest term (optimal truncation). F, G are Si's f, g without the (-1)^k. The odd-symmetry reduction Shi(-z) = -Shi(z) folds negative real / left-half-plane inputs onto Re >= 0, keeping the asymptotic within its valid sector. On the real axis the bare form above is exact (real, no constant); off it, matching Mathematica's Series[SinhIntegral[z],{z,Infinity,k}], a Stokes constant restores the analytic value:

```text
  Shi(z) = B(z) + i (Pi/2) sign(Im z),   B(z) = cosh(z) F(z) + sinh(z) G(z).
```

(B is odd; on the imaginary axis B ~ O(1/z) vanishes and the i Pi/2 sign(Im z) is all that survives, giving Shi(+- i Infinity) = +- i Pi/2.)

Machine-real results that overflow a C double (e.g. SinhIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[SinIntegral](../../special-functions/SinIntegral/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_sinhintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sinhintegral.c)
