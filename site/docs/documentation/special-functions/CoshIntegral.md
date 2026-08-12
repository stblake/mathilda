# CoshIntegral

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CoshIntegral[z]`**

gives the hyperbolic cosine integral Chi(z) = EulerGamma + Log\[z\] + Integral\_0^z (Cosh\[t\] - 1)/t dt.

**`CoshIntegral[0] = -Infinity, CoshIntegral[Infinity] = Infinity,`**

**`CoshIntegral[+-I Infinity] = +-I Pi/2.`**

<details>
<summary>Notes</summary>

Has a logarithmic singularity at 0 and a branch cut on (-Infinity, 0\]. Real and complex inputs evaluate numerically at machine or arbitrary (MPFR) precision; D\[CoshIntegral\[z\], z\] = Cosh\[z\]/z. Listable.

</details>

## Examples

_No verified examples yet for this function._

## Algorithm

Mathilda -- the hyperbolic cosine integral

```text
  Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.

  CoshIntegral[z]   Chi(z)
```

Chi has a logarithmic singularity at z = 0 and a branch cut running along the

```text
negative real axis (-Infinity, 0].  It is the imaginary-axis sibling of the
cosine integral: Chi(z) = Ci(i z) - i Pi/2.  Evaluation is layered so each
```

argument takes the cheapest route:

```text
  exact special values     ->  -Infinity, Infinity, +-I Pi/2, Indeterminate
  machine real             ->  MPFR series/asymptotic at 53 bits
  arbitrary real           ->  the same, at the input precision
  complex (any precision)  ->  the ncpx series/asymptotic with guard bits
  everything else          ->  stays symbolic (return NULL)
```

The convergent Maclaurin series (the trig series with the alternating sign removed) is

```text
  Chi(z) = EulerGamma + Log(z) + Sum_{k>=1} z^(2k) / (2k (2k)!),
```

valid on the principal branch for all z != 0 -- the principal Log(z) supplies the correct +-i Pi jump across the cut, so no folding is needed for the

```text
convergent path.  Every series term is positive, so on the real axis there is
NO catastrophic cancellation and only a small fixed guard is needed.  For large
|z| the convergent series is infeasible; there we use the asymptotic expansion

  Chi(z) ~ sinh(z) F(z) + cosh(z) G(z)   (+ Stokes constant, see below),
    F(z) ~ Sum (2k)!   / z^(2k+1) = 1/z + 2!/z^3 + ...,
    G(z) ~ Sum (2k+1)! / z^(2k+2) = 1/z^2 + 3!/z^4 + ...,
```

(the same F, G SinhIntegral computes; only the combination differs -- Shi uses

```text
cosh F + sinh G).  The bare part sinh F + cosh G is even and real on the
positive real axis (no constant).  Off the real axis, matching Mathematica's
```

Series[CoshIntegral[z],{z,Infinity,k}], a piecewise Stokes term restores the principal branch:

```text
  Chi(z) = B2(z) + i K,   B2(z) = sinh(z) F(z) + cosh(z) G(z),
    K = -Pi/2                       for Im(z) < 0,
    K = Pi sgn+(Re z) - Pi/2        for Im(z) > 0   (sgn+(0) = +1, from above).
```

(On the imaginary axis Im(z) > 0, Re(z) = 0 gives K = +Pi/2, so

```text
Chi(+- i Infinity) = +- i Pi/2.)  A negative real x is handled by the real
```

path as the from-above branch value Complex[Chi(|x|), Pi], matching the log jump; and Chi(-Infinity) -> Infinity (the real part dominates).

Machine-real results that overflow a C double (e.g. CoshIntegral[10.^6] ~ 1.5*10^434288) are emitted as a 53-bit MPFR real, whose exponent range easily holds them, rather than as an inf-valued double.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[CosIntegral](../../special-functions/CosIntegral/), [Log](../../elementary-functions/Log/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_coshintegral.c`](https://github.com/stblake/mathilda/blob/main/tests/test_coshintegral.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)
- Tests: [`tests/test_series.c`](https://github.com/stblake/mathilda/blob/main/tests/test_series.c)
