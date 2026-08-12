# Gamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Gamma[z]`**

is the Euler gamma function Gamma(z).

**`Gamma[a, z]`**

is the upper incomplete gamma function Gamma(a, z).

**`Gamma[a, z0, z1]`**

is the generalized incomplete gamma Gamma(a, z0) - Gamma(a, z1).

<details>
<summary>Notes</summary>

Integer and half-integer arguments reduce to exact values ((z-1)!, and rational multiples of Sqrt\[Pi\]); non-positive integers give ComplexInfinity. Machine and arbitrary-precision (MPFR) real inputs evaluate numerically, as do machine-precision complex inputs. Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= Gamma[5]
Out[1]= 24
```

```mathematica
In[1]:= Gamma[7/2]
Out[1]= 15/8 Sqrt[Pi]
```

```mathematica
In[1]:= Gamma[-1/2]
Out[1]= -2 Sqrt[Pi]
```

```mathematica
In[1]:= N[Gamma[1/3], 40]
Out[1]= 2.6789385347077476336556929409746776441289
```

```mathematica
In[1]:= N[Gamma[3 + 4 I], 20]
Out[1]= 0.00522553847136921419473 - 0.172547079294300187719*I
```

## Algorithm

Mathilda -- the Gamma function family.

```text
  Gamma[z]          Euler gamma function   Gamma(z) = Int_0^Inf t^(z-1) e^-t dt
  Gamma[a, z]       upper incomplete gamma Gamma(a,z) = Int_z^Inf t^(a-1) e^-t dt
  Gamma[a, z0, z1]  generalized incomplete = Gamma(a,z0) - Gamma(a,z1)
```

Evaluation is layered so each kind of argument takes the cheapest exact or fastest numeric route available:

```text
  exact integer / half-integer  ->  (z-1)! via the Factorial machinery
                                     (exact, BigInt, or rational*Sqrt[Pi])
  machine real        ->  libm   tgamma
  machine complex     ->  Lanczos approximation (double complex)
  arbitrary real      ->  MPFR   mpfr_gamma / mpfr_gamma_inc
  everything else     ->  stays symbolic (return NULL)
```

Arbitrary-precision *complex* gamma is deliberately left symbolic: a fixed-coefficient Lanczos series only carries ~15 correct digits, so emitting it as an MPFR value would advertise a precision it does not have. Reporting the input unevaluated is the honest behaviour.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

- Exact reductions for `Gamma[z]`:
  - Positive integers: `Gamma[n] = (n-1)!` (exact, with GMP BigInt for large `n`).
  - Non-positive integers are poles: `Gamma[0]`, `Gamma[-n]` → `ComplexInfinity`.
  - Half-integers reduce to rational multiples of `Sqrt[Pi]`, e.g.
    `Gamma[1/2] = Sqrt[Pi]`, `Gamma[5/2] = 3/4 Sqrt[Pi]`,
    `Gamma[-1/2] = -2 Sqrt[Pi]` (via the Factorial functional equation).
  - `Gamma[Infinity]` → `Infinity`, `Gamma[-Infinity]` → `Indeterminate`,
    `Gamma[ComplexInfinity]` → `ComplexInfinity`.
- Exact reductions for the incomplete form:
  - `Gamma[a, 0] = Gamma[a]`, `Gamma[a, Infinity] = 0`.
  - **Positive integer first argument** reduces to its finite closed form
    `Gamma[n, z] = (n-1)! e^-z Σ_{k=0}^{n-1} z^k/k!` for symbolic or exact `z`,
    e.g. `Gamma[1, z] = E^-z`, `Gamma[2, x] = (1 + x) E^-x`,
    `Gamma[3, x] = (2 + 2 x + x^2) E^-x`, and exact `Gamma[2, 3] = 4/E^3`.
- Numeric evaluation:
  - Machine-precision real → libm `tgamma`; machine-precision complex →
    Lanczos approximation, e.g. `Gamma[2.3 + I] = 0.719141 + 0.540614 I`.
  - Arbitrary precision (MPFR) real → `mpfr_gamma`, output precision tracking
    the input, e.g. `N[Gamma[22/10], 50]` and `Gamma[2.2`200]`.
  - **Arbitrary-precision complex** `Gamma[z]` → Spouge's approximation,
    whose coefficients are computed at runtime to the requested precision
    (reflection for `Re(z) < 1/2`), e.g.
    `N[Gamma[I], 50] = -0.15494982830181068512… − 0.49801566811835604271… I`.
  - Incomplete real (machine or MPFR) → `mpfr_gamma_inc`, e.g.
    `Gamma[1.5, 7.5] = 0.00160996`, `Gamma[1, 1.1, 2.2] = 0.222068`.
  - **Incomplete complex** `Gamma[a, z]` (machine or arbitrary precision) →
    a lower-incomplete series (`Re(z) < Re(a)+1`) or a Lentz continued
    fraction otherwise, e.g. `Gamma[2.0, 1 + I] = (2 + I) e^-(1+I)` and
    `N[Gamma[3/2, 2 + I], 30] = 0.160487401929263240… − 0.176588715957602346… I`.
- Derivatives: `D[Gamma[a, z], z] = -z^(a-1) E^-z` (chain rule on both
  arguments; the `a`-derivative is the generic `Derivative[1,0][Gamma][a,z]`).
  `D[Gamma[z], z] = Gamma[z] PolyGamma[0, z]`, so higher derivatives compose
  through `PolyGamma`, e.g.
  `D[Gamma[z], {z, 2}] = Gamma[z] PolyGamma[1, z] + Gamma[z] PolyGamma[0, z]^2`.
- All other arguments (e.g. `Gamma[1/3]`, `Gamma[x]`, `Gamma[a, z]`,
  exact non-integer `Gamma[3/2, z]`, exact complex `Gamma[3/2, I]`) stay
  unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[PolyGamma](../../special-functions/PolyGamma/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_beta.c`](https://github.com/stblake/mathilda/blob/main/tests/test_beta.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`Gamma` is the Euler gamma function, the analytic continuation of the
factorial: `Gamma[n] = (n-1)!`, so `Gamma[5] = 24`. Half-integer arguments
collapse to exact rational multiples of `Sqrt[Pi]` — `Gamma[7/2] = 15/8 Sqrt[Pi]`
— and this continues through the poles at non-positive integers into negative
half-integers, where `Gamma[-1/2] = -2 Sqrt[Pi]`. For arguments with no
closed form it evaluates to arbitrary precision via MPFR (`Gamma[1/3]` to 40
digits) and across the complex plane (`Gamma[3 + 4 I]`). The two- and
three-argument forms give the upper incomplete and generalized incomplete gamma
functions.
