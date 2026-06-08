# Rationalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Rationalize[x]
    converts an approximate number x to a nearby rational with small denominator.
Rationalize[x, dx]
    yields the rational number with smallest denominator that lies within dx of x.

Rationalize[x] yields x unchanged if there is no rational number close enough to x to satisfy |p/q - x| < c/q^2, with c = 10^-4.
Rationalize[x, dx] works with exact numbers x: the value is first numericalised, then rationalised.
Rationalize[x, 0] forces conversion of any inexact number x to rational form, using a tolerance derived from the precision of x.
Rationalize threads over compound expressions and Complex[re, im], so e.g. Rationalize[1.2 + 6.7 x] gives 6/5 + (67 x)/10.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_rationalize` dispatches on arity into `internal_rationalize_expr` with one of three modes. `Rationalize[x]` (`RATIONALIZE_DEFAULT`) walks the convergents of the continued-fraction expansion of `x` and returns the first `p/q` whose error satisfies `|p/q − x| < 10^-4 / q^2` (Mathematica's threshold); if none qualify it returns the value unchanged. `Rationalize[x, dx]` with `dx > 0` (`RATIONALIZE_TOLERANCE`) solves the "simplest rational in an interval" problem by the Stern–Brocot / continued-fraction descent: it extracts integer parts of the lower and upper bounds until they diverge, at which point the smallest integer in the lower interval gives the minimal-denominator `p/q` with `|p/q − x| ≤ dx`. `dx == 0` (`RATIONALIZE_ZERO`) is treated as a thin ulp-derived tolerance so e.g. `Rationalize[N[Pi], 0]` yields `245850922/78256779` rather than the bit-exact dyadic. GMP is used throughout so denominators exceeding `int64_t` are handled. A negative or non-numeric `dx` returns `NULL` (left unevaluated).

**Threading.** Rationalize descends through expression structure: a compound subexpression that is itself `NumericQ` (e.g. `Sqrt[2]`, `Exp[Sqrt[2]]`) is numericalised via the `N` pipeline and then rationalised end-to-end; otherwise the head and arguments are recursively rationalised and re-evaluated so `Rationalize[1.2 + 6.7 x]` becomes `6/5 + (67 x)/10`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/rationalize.c`](https://github.com/stblake/mathilda/blob/main/src/rationalize.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Rationalize[0.5]
Out[1]= 1/2

In[2]:= Rationalize[N[Pi], 10^-4]
Out[2]= 333/106

In[3]:= Rationalize[1.2 + 6.7 x]
Out[3]= 6/5 + 67/10 x
```

### Notes

`Rationalize[x]` finds a nearby rational with small denominator (within `c/q^2`, `c = 10^-4`); the two-argument form `Rationalize[x, dx]` returns the smallest-denominator rational within `dx` of `x`. It threads over compound expressions.
