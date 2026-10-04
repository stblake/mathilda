# InverseErf

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InverseErf[s]`**

gives the inverse error function: the z solving s = Erf\[z\].

**`InverseErf[z0, s]`**

gives the inverse of the generalized error function Erf\[z0, z\].

**`InverseErf[0] = 0, InverseErf[1] = Infinity, InverseErf[-1] = -Infinity.`**

<details>
<summary>Notes</summary>

Odd in s. Numerical values are given only for real s in \[-1, 1\], at machine or arbitrary (MPFR) precision; D\[InverseErf\[z\], z\] = (Sqrt\[Pi\]/2) E^(InverseErf\[z\]^2). Listable.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= InverseErf[0]
Out[1]= 0

In[2]:= N[InverseErf[1/2], 40]
Out[2]= 0.47693627620446987338141835364313055980899

In[3]:= Series[InverseErf[x], {x, 0, 7}]
Out[3]= 1/2 Sqrt[Pi] x + 1/24 Pi^(3/2) x^3 + 7/960 Pi^(5/2) x^5 + 127/80640 Pi^(7/2) x^7 + O[x]^8

In[4]:= D[InverseErf[z], z]
Out[4]= 1/2 Sqrt[Pi] E^InverseErf[z]^2

In[5]:= N[Sqrt[2] InverseErf[2 (95/100) - 1], 30]
Out[5]= 1.644853626951472714863848907992
```

## Algorithm

Mathilda -- the inverse error function.

```text
  InverseErf[s]       inverse of erf: the z solving s = erf(z)
  InverseErf[z0, s]   inverse of the generalized error function:
                      the z solving s = Erf[z0, z] = erf(z) - erf(z0)
```

Per Mathematica, explicit numerical values are produced only for *real* s in [-1, 1]; complex and out-of-domain (|s| > 1) inputs stay symbolic. Evaluation is therefore much simpler than Erf -- no complex machinery:

```text
  exact special values   ->  InverseErf[0]=0, InverseErf[1]=Infinity,
                             InverseErf[-1]=-Infinity
  symbolic odd argument  ->  InverseErf[-x] = -InverseErf[x]
  machine real (|s|<1)   ->  Winitzki seed + Newton polish on libm erf
  arbitrary real (|s|<1) ->  Newton with precision doubling on mpfr_erf
  everything else        ->  stays symbolic (return NULL)
```

Neither C99 libm nor MPFR ships an inverse-erf, so the kernel is Newton's method on f(z) = erf(z) - s, with f'(z) = (2/sqrt(pi)) e^{-z^2}, i.e.

```text
  z <- z - (erf(z) - s) (sqrt(pi)/2) e^{z^2}.
```

Newton converges quadratically; the MPFR path doubles the working precision each step so the final full-precision erf dominates the cost.

The derivative D[InverseErf[z], z] = (sqrt(pi)/2) e^{InverseErf[z]^2} lives in src/calculus/deriv.c; Series follows from it via the generic Taylor-via-D fallback.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_inverf` solves `s = erf(z)`. Per Mathematica, explicit
values are produced only for *real* `s` in `[-1, 1]`; complex and `|s| > 1`
inputs stay symbolic. Exact: `0 -> 0`, `1 -> Infinity`, `-1 -> -Infinity`; a
symbolic negative-leading `Times` folds by odd symmetry `InverseErf[-x] ->
-InverseErf[x]`. Machine real `|s| < 1` uses `inverf_double`: a Winitzki
closed-form seed polished by four Newton steps on the libm `erf`,
`z <- z - (erf(z) - s)(sqrt(pi)/2) e^{z^2}`. Arbitrary-precision real uses Newton
with precision doubling on `mpfr_erf`, the final full-precision `erf` dominating
cost. The two-argument `InverseErf[z0, s]` rewrites to
`InverseErf[s + Erf[z0]]`, staying symbolic if `Erf[z0]` does not reduce. The
derivative `(sqrt(pi)/2) e^{InverseErf[z]^2}` lives in `calculus/deriv.c`; Series
follows via Taylor-via-`D`.

**Data structures.** `Expr`; `double` (`inverf_double`, exposed for ND and the
`Compile[]` VM); `mpfr_t` for the arbitrary-precision Newton loop. ND: real-only
unary kernel `NDKU_InverseErf = { NULL, ndk_InverseErf_r, ... }` (no complex
kernel), registered `REG_U`, so `packed_aware`. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Newton converges quadratically; the MPFR cost tracks
the final full-precision `erf`. Domain is the real interval `[-1, 1]`; outside
it (and for complex or exact-rational `s`) the call stays symbolic.
`Compile[]` lowers at both scalar and rank-1 array shapes (`Compiled -> True`).

- Exact special values: `InverseErf[0] = 0`, `InverseErf[1] = Infinity`,
  `InverseErf[-1] = -Infinity` (likewise `InverseErf[1.] = Infinity`), plus
  `Indeterminate` passes through.
- Odd symmetry for symbolic arguments: `InverseErf[-x] = -InverseErf[x]`.
- Numeric evaluation (neither C99 libm nor MPFR ships an inverse-erf, so the
  kernel is Newton's iteration on `f(z) = erf(z) − s`,
  `z ← z − (erf(z) − s)(Sqrt[Pi]/2) e^{z²}`):
  - Machine-precision real → a Winitzki seed polished by Newton on libm `erf`,
    e.g. `InverseErf[0.6] = 0.595116`, `InverseErf[1/{2.,3.,4.,5.}] =
    {0.476936, 0.30457, 0.225312, 0.179143}`.
  - Arbitrary precision (MPFR) real → Newton with precision doubling on
    `mpfr_erf`, output precision tracking the input, e.g.
    `N[InverseErf[33/100], 50] = 0.30133214613370582612850271815839477396582428282853`.
- Two-argument form: `InverseErf[0.4, 0.2] = 0.631776`; if `Erf[z0]` does not
  reduce the call stays in two-argument form, while a reducible `z0` collapses to
  the one-argument inverse, e.g. `InverseErf[0, 1.3] = InverseErf[1.3]`.
- Derivative: `D[InverseErf[z], z] = (Sqrt[Pi]/2) E^(InverseErf[z]²)` (chain rule
  applies), so the origin Taylor series follows from the generic `D`-based
  fallback, e.g. `Series[InverseErf[x], {x, 0, 8}] = (Sqrt[Pi]/2) x +
  (Pi^(3/2)/24) x³ + (7 Pi^(5/2)/960) x⁵ + (127 Pi^(7/2)/80640) x⁷ + O[x]^9`.
- All other arguments (symbolic `InverseErf[x]`, exact `InverseErf[1/2]`,
  out-of-domain `InverseErf[2]`) stay unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [D](../../calculus/D/)

- DLMF §7.17 — inverse error functions.
- S. Winitzki, A handy approximation for the error function and its inverse (2008) — the closed-form seed.
- Source: [`src/special_functions/inverf.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/inverf.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_inverf.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inverf.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`InverseErf[s]` returns the `z` solving `Erf[z] == s`. It is odd in `s`, with
`InverseErf[0] = 0`, `InverseErf[1] = Infinity`, `InverseErf[-1] = -Infinity`.
Numerical values are produced only for real `s` in `[-1, 1]`, at machine or
arbitrary (MPFR) precision; symbolic arguments are returned unevaluated.
