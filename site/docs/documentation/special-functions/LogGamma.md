# LogGamma

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LogGamma[z]`**

gives the log-gamma function log(Gamma(z)), analytic except for a branch

<details>
<summary>Notes</summary>

cut on the negative reals. Exact at integer and half-integer z (with the negative-axis branch term), divergent (Infinity) at non-positive integers, and evaluated numerically for real and complex z at machine or arbitrary (MPFR) precision. D\[LogGamma\[z\], z\] is PolyGamma\[0, z\]. Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= LogGamma[5]
Out[1]= Log[24]

In[2]:= LogGamma[1/2]
Out[2]= Log[Sqrt[Pi]]

In[3]:= D[LogGamma[z], z]
Out[3]= PolyGamma[0, z]

In[4]:= N[LogGamma[100], 40]
Out[4]= 359.13420536957539877604401046028690961264

In[5]:= N[LogGamma[1 + I], 30]
Out[5]= -0.6509231993018563388852168315042 - 0.3016403204675331978875316577968*I
```

## Algorithm

Mathilda -- LogGamma[z], the log-gamma function log(Gamma(z)).

LogGamma is analytic throughout the complex plane except for a branch cut on the negative real axis. It is the analytic continuation of log(Gamma(z)) and is *not* Log[Gamma[z]] (which has a more involved branch structure). Because it never overflows where Gamma does, it is the natural primitive for factorial ratios and asymptotics.

Evaluation is layered, mirroring Gamma (gamma.c), so each kind of argument takes the cheapest exact or fastest numeric route:

```text
  exact integer n >= 1          ->  Log[(n-1)!]   (exact, via Gamma's machinery)
  exact half-integer            ->  Log of the exact Sqrt[Pi] form, with the
                                    branch term -Ceiling[-z] Pi I for z < 0
  non-positive integer          ->  Infinity      (pole)
  symbolic infinities           ->  Infinity / Indeterminate / ComplexInfinity
  machine real z > 0            ->  libm   lgamma
  machine real z < 0 (non-int)  ->  lgamma + branch term (complex result)
  arbitrary real                ->  MPFR   mpfr_lgamma  (+ branch term)
  machine complex               ->  Lanczos log-gamma (double complex)
  arbitrary complex             ->  Spouge log-gamma   (MPFR, runtime coeffs)
  everything else               ->  stays symbolic (return NULL)
```

The imaginary branch term on the negative real axis is Im = -Pi Ceiling[-z], taken from above (Mathematica's convention); this is exact for the symbolic half-integer path and is reused by the numeric real paths.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

- **Exact closed forms.** Integers reduce as `LogGamma[n] = Log[(n-1)!]`
  (`LogGamma[5] = Log[24]`); positive half-integers give `Log` of the exact
  `Sqrt[Pi]` form; negative half-integers carry the branch term
  `-Ceiling[-z] Pi I`, e.g. `LogGamma[-3/2] = -2 I Pi + Log[(4 Sqrt[Pi])/3]`.
- **Poles.** Non-positive integers diverge: `LogGamma[0] = LogGamma[-1] = … =
  Infinity`.
- **Symbolic infinities.** `LogGamma[Infinity] = Infinity`,
  `LogGamma[-Infinity] = Indeterminate`, `LogGamma[I Infinity] =
  LogGamma[ComplexInfinity] = ComplexInfinity`.
- **Branch continuity across the reflection.** For `Re[z] < 1/2` both numeric
  paths use `LogGamma[z] = Log[Pi] - Log[Sin[Pi z]] - LogGamma[1-z]`, where the
  log of the sine must be the *continued* one. Factoring
  `Sin[Pi z] = E^(-I Pi z) (E^(2 I Pi z) - 1)/(2 I)` puts the whole winding in
  the exact `-I Pi z` term and leaves a factor near `-1` whose principal log is
  safe. At `Im[z] = 0` this is the limit from above, matching the
  `Im = -Ceiling[-z] Pi` convention used on the real axis. Before 2026-07-27 a
  principal `Log[Pi/Sin[Pi z]]` was used instead: correct only in the strip
  `-1 < Re[z] < 0`, and short by a multiple of `2 Pi I` outside it, so
  `LogGamma[-4.5 + 3. I]` equalled `Log[Gamma[-4.5 + 3. I]]` rather than the
  continuation.
- **Numerics.** Machine real via `lgamma`; arbitrary-precision real via MPFR
  `lgamma`; machine complex via a Lanczos log-series; arbitrary-precision
  complex via a Stirling series with argument reduction (the continuous
  branch — its imaginary part grows past π where `Log[Gamma]` would wrap).
  Negative real arguments return the complex value `log|Γ(z)| - Pi Ceiling[-z] I`.
- **Derivative.** `D[LogGamma[z], z] = PolyGamma[0, z]`; higher derivatives
  raise the `PolyGamma` order. Produced by `PolyGamma[-1, z]`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Gamma](../../special-functions/Gamma/), [Log](../../elementary-functions/Log/), [PolyGamma](../../special-functions/PolyGamma/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_fullsimplify.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fullsimplify.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)

## Notes & additional examples

### Notes

`LogGamma[z]` is `log(Gamma(z))`, analytic except for a branch cut on the negative reals. It is exact at integer and half-integer arguments (`LogGamma[5]` is `Log[4!] = Log[24]`), divergent at non-positive integers, and evaluates numerically for real or complex `z` at machine or arbitrary (MPFR) precision. Its derivative is `PolyGamma[0, z]`. Unlike `Log[Gamma[z]]`, `LogGamma` tracks the correct sheet, which matters for large or complex arguments where `Gamma` overflows. Listable.
