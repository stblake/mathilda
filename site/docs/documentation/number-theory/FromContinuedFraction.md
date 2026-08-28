# FromContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromContinuedFraction[{a1, a2, ..., an}] reconstructs a1 + 1/(a2 + 1/(a3 + ... + 1/an)), the convergent of a continued fraction, in nested (un-expanded) form; the inverse of ContinuedFraction.`**

**`FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]`**

gives the exact quadratic irrational whose continued-fraction terms begin with the ai then cycle through the bi forever; all ai and bi must be integers. FromContinuedFraction\[{}\] is 0.  The terms may be symbolic.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FromContinuedFraction[{2, 1, 3, 4}]
Out[1]= 47/17

In[2]:= FromContinuedFraction[{a, b, c, d}]
Out[2]= (1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)

In[3]:= FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]
Out[3]= Sqrt[71]

In[4]:= FromContinuedFraction[{{1, 2, 3, 4}}]
Out[4]= 1/15 (9 + 2 Sqrt[39])

In[5]:= FromContinuedFraction[ContinuedFraction[Pi, 3]]
Out[5]= 333/106
```

### Applications (4)

```mathematica
In[6]:= FromContinuedFraction[{3, 7, 15, 1, 292}]
Out[6]= 103993/33102

In[7]:= FromContinuedFraction[{1, {1}}]
Out[7]= 1/2 (1 + Sqrt[5])

In[8]:= FromContinuedFraction[{1, 2, {2}}]
Out[8]= Sqrt[2]

In[9]:= FromContinuedFraction[{a, b, c}]
Out[9]= (a + (1 + a b) c)/(1 + b c)
```

## Algorithm

contfrac.c — ContinuedFraction[x] and ContinuedFraction[x, n].

Computes the *simple* continued-fraction expansion of x:

```text
    {a1, a2, a3, ...}  <->  a1 + 1/(a2 + 1/(a3 + ...))
```

Four input regimes are handled, each with the appropriate exact or precision-aware algorithm:

1. Exact rationals (Integer, BigInt, Rational[p, q]).

```text
    The Euclidean algorithm with floor division produces the canonical
    terminating expansion (last term >= 2, never the {..., k-1, 1} form).
    A finite rational yields finitely many terms, so ContinuedFraction[x, n]
    may return fewer than n terms.
```

2. Quadratic irrationals of the form Sqrt[D], D a positive non-square

```text
    integer.  The classic periodic surd recurrence
        m_{i+1} = a_i q_i - m_i,
        q_{i+1} = (D - m_{i+1}^2) / q_i,
        a_{i+1} = floor((a0 + m_{i+1}) / q_{i+1})
    is purely periodic after a0 for an integer square root, so the period
    is detected when the state (m, q) first repeats.  Without n the result
    is {a0, {b1, ..., bk}} (the bracketed list is the repeating block);
    with n the periodic sequence is unrolled to exactly n terms.
    General quadratic irrationals (e.g. (1 + Sqrt[5])/2, Sqrt[2/3]) are not
    recognised symbolically here — supply an explicit n to use the numeric
    path instead.

 3. Inexact reals (machine Real, or arbitrary-precision MPFR).  Terms are
    extracted by repeated reciprocation while tracking the absolute
    uncertainty of the running value (initially |x| * 2^-prec).  Expansion
    stops as soon as the integer part of the value is no longer determined
    by the available precision — i.e. ContinuedFraction stops when it runs
    out of precision.  The iteration carries 64 guard bits so arithmetic
    round-off stays far below the modelled input uncertainty.
```

4. Exact symbolic reals with an explicit n (Pi, Sqrt[E], Exp[Pi Sqrt[163]],

```text
    ...).  The value is numericised via N[x, digits] at successively
    doubled precision until n reliable terms are obtained.
```

Anything else (an exact non-rational, non-quadratic with no n, or a non-real numeric value) is left unevaluated (the builtin returns NULL).

Attributes: Listable, Protected.

## Implementation notes

**Algorithm.** `builtin_from_continued_fraction` reconstructs the value from a `List` of terms. A non-periodic list uses the standard convergent recurrence `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}` (`fcf_simple`), evaluating each step so numeric terms collapse and symbolic terms stay in convergent form; the result is `h_{n-1}/k_{n-1}`. A trailing sub-list marks the cyclic (period) block, requiring all integer terms: `fcf_periodic` builds the period's convergents in GMP, forms the quadratic `A x^2 + B x + C = 0` for the purely-periodic tail (with `A=k_{k-1}`, `B=k_{k-2}-h_{k-1}`, `C=-h_{k-2}`), solves via the discriminant with largest-square extraction (`fcf_extract_square`, trial division up to `10^6`), then applies the leading terms as a Möbius transform `(Hx+H')/(Kx+K')` and rationalises to `(P + Q Sqrt[R])/S` (`fcf_qirr_to_expr`).

**Data structures.** Convergents are GMP `mpz_t` registers; symbolic non-periodic convergents are `Expr` trees via `eval_and_free`. Output is an Integer/Rational or a rationalised quadratic-irrational expression.

**Complexity / limits.** Linear in the number of terms. The periodic path requires exact integer terms; a residual `p^2 q` with both `p, q > 10^6` is left un-reduced (astronomically unlikely for reconstructed CF data).

- `Protected` (not `Listable` — the argument is the whole term list).
- The `ai` of the finite form may be **symbolic**; the result is the convergent
  `h_n / k_n` built from the fundamental recurrence
  `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}`, kept in nested
  (un-expanded) form — `Together` collapses it to a flat rational.
- The **periodic** form requires all `ai` and `bi` to be integers. The purely
  periodic tail solves the quadratic
  `k_{k-1} x^2 + (k_{k-2} - h_{k-1}) x - h_{k-2} = 0` (h, k the period's
  convergents); its positive root is then pushed through the leading terms by a
  Möbius transform and rationalised to a single `(P + Q Sqrt[R]) / S` in lowest
  terms.
- `FromContinuedFraction[{}]` is `0`; `FromContinuedFraction[{x}]` is `x`.
- Left unevaluated for a non-list argument, a sub-list anywhere but last, an
  empty period block, or non-integer terms in a periodic form.

**Attributes:** `Protected`.

## References

**See also:** [ContinuedFraction](../../number-theory/ContinuedFraction/), [Together](../../algebra/Together/)

- A. Ya. Khinchin, *Continued Fractions*, Dover, 1997 — convergents and their approximation properties.
- Source: [`src/contfrac.c`](https://github.com/stblake/mathilda/blob/main/src/contfrac.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_contfrac.c`](https://github.com/stblake/mathilda/blob/main/tests/test_contfrac.c)

## Notes & additional examples

### Folding a continued fraction back

`FromContinuedFraction` is the inverse of [`ContinuedFraction`](ContinuedFraction.md): it
folds a list of partial quotients into the rational (or, from a periodic block, the exact
quadratic irrational) it represents. Truncating an expansion and folding it back yields a
*convergent* — the best rational approximation of its size — as with `{3, 7, 15, 1, 292}`
giving `103993/33102 ≈ π`.

### Notes

For a finite list, `FromContinuedFraction` reconstructs the rational (or
symbolic) convergent `a1 + 1/(a2 + 1/(a3 + ...))`. The classic terms
`{3, 7, 15, 1, 292}` give `103993/33102`, the celebrated convergent of Pi.

A trailing sublist marks the periodic part of a *quadratic irrational*: the
purely periodic `{1, {1}}` returns the golden ratio `(1 + Sqrt[5])/2`, and
`{1, 2, {2}}` recovers `Sqrt[2]` from its eventually-periodic expansion
`[1; 2, 2, 2, ...]`. With symbolic terms the result is the exact nested form,
left un-expanded. `FromContinuedFraction` is the inverse of `ContinuedFraction`.
