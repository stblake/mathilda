# ContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ContinuedFraction[x, n] gives the first n terms of the continued-fraction expansion of x, the list {a1, a2, ...} standing for a1 + 1/(a2 + 1/(a3 + ...)); truncating it gives the convergents, the best rational approximations to x.`**

**`ContinuedFraction[x]`**

gives all terms determinable from the precision of x.

<details>
<summary>Notes</summary>

Exact rationals give a finite (canonical, last term \>= 2) expansion.  For Sqrt\[d\] with d a non-square integer the no-count form returns {a1, ..., {b1, ...}}, the bracketed block repeating cyclically -- every quadratic irrational is eventually periodic (Lagrange).  Inexact Real / MPFR inputs yield terms only as far as the precision determines them.  ContinuedFraction is Listable.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ContinuedFraction[47/17]
Out[1]= {2, 1, 3, 4}

In[2]:= ContinuedFraction[Sqrt[13]]
Out[2]= {3, {1, 1, 1, 1, 6}}

In[3]:= ContinuedFraction[Pi, 20]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}

In[4]:= ContinuedFraction[N[Pi]]
Out[4]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}

In[5]:= ContinuedFraction[Exp[Pi Sqrt[163]], 10]
Out[5]= {262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}
```

### Applications (5)

```mathematica
In[6]:= ContinuedFraction[123/47]
Out[6]= {2, 1, 1, 1, 1, 1, 1, 3}

In[7]:= ContinuedFraction[Sqrt[2]]
Out[7]= {1, {2}}

In[8]:= ContinuedFraction[Sqrt[7]]
Out[8]= {2, {1, 1, 1, 4}}

In[9]:= ContinuedFraction[N[Pi, 40], 12]
Out[9]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1}

In[10]:= FromContinuedFraction[{1, 2, 2, 2, 2, 2}]
Out[10]= 99/70
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

**Algorithm.** `builtin_continued_fraction` computes the *simple* continued-fraction expansion, dispatching on input regime: (1) **exact rationals** use the Euclidean algorithm with floor division (`cf_rational`), producing the canonical terminating form (last term `>= 2`); (2) **`Sqrt[D]`** for a non-square positive integer `D` uses the classic periodic-surd recurrence `m'=a q-m, q'=(D-m'^2)/q, a'=floor((a0+m')/q')`, detecting the period when the `(m, q)` state first repeats (`cf_sqrt_period`), and returns `{a0, {period...}}` (or unrolls to `n` terms); (3) **inexact reals** (machine or MPFR) extract terms by repeated reciprocation while tracking absolute uncertainty (`|x| 2^-prec` plus 64 guard bits), stopping when the integer part is no longer determined by available precision (`cf_inexact`); (4) **exact symbolic reals with explicit `n`** (Pi, etc.) are numericised via `N[x, digits]` at doubling precision until `n` terms stabilise (`cf_exact_numeric`).

**Data structures.** Terms accumulate in a growable `TermVec` of GMP `mpz_t`; the MPFR path uses `mpfr_t` working registers; output is a `List` of integers (with a nested `List` for the repeating block in the unbounded surd case).

**Complexity / limits.** Rationals terminate in `O(log)` Euclidean steps. General quadratic irrationals beyond bare `Sqrt[D]` (e.g. `(1+Sqrt[5])/2`) are not recognised symbolically — supply an explicit `n` to use the numeric path. Inexact inputs stop at the precision floor.

- `Protected`, `Listable`.
- **Exact rationals** (Integer / BigInt / Rational) use the Euclidean
  algorithm and return the canonical terminating form (last term `>= 2`,
  never `{..., k-1, 1}`). A finite rational may yield fewer than `n` terms.
- **Quadratic surds** `Sqrt[d]` with `d` a non-square integer use the
  periodic surd recurrence. Without a count the result is
  `{a1, ..., {b1, ...}}`, where the bracketed block repeats cyclically; with
  a count the periodic sequence is unrolled to exactly `n` terms. (General
  quadratic irrationals are not recognised symbolically — pass an explicit
  `n` to use the numeric path.) The no-count form is declined for a `d` whose
  period would be astronomically long.
- **Inexact reals** (machine `Real` or arbitrary-precision MPFR) yield terms
  only as far as the input precision determines them, tracking the value's
  uncertainty and stopping when the next term is no longer pinned down.
- **Exact symbolic reals with an explicit `n`** (e.g. `Pi`, `Sqrt[E]`,
  `Exp[Pi Sqrt[163]]`) are numericised at adaptively increasing precision
  until `n` terms are confirmed by two consecutive evaluations.
- Left unevaluated for an exact non-rational, non-quadratic value with no
  count, or a non-real numeric value.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Pi](../../mathematical-constants/Pi/)

- A. Ya. Khinchin, *Continued Fractions*, Dover, 1997 — the classical theory of convergents and best approximations.
- G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 6th ed., Oxford University Press, 2008 — Chapter X, continued fractions.
- Source: [`src/contfrac.c`](https://github.com/stblake/mathilda/blob/main/src/contfrac.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_contfrac.c`](https://github.com/stblake/mathilda/blob/main/tests/test_contfrac.c)

## Notes & additional examples

### Convergents and periodicity

The partial quotients of `ContinuedFraction[p/q]` are exactly the quotients of the Euclidean
algorithm on `p` and `q`, so a rational has a *finite* expansion. Truncating the expansion
and folding it back (`FromContinuedFraction`) gives the *convergents* — the best rational
approximations to `x` for their size, which is why `{3, 7, 15, 1}` recovers `355/113 ≈ π`. A
quadratic irrational such as `Sqrt[7]` has an *eventually periodic* expansion (Lagrange's
theorem), returned in the no-count form as a bracketed repeating block.

### Notes

For an exact rational the expansion is finite and canonical (last term `>= 2`). For `Sqrt[d]` with `d` a non-square integer the no-count form returns the periodic block in braces, `{a1, {b1, ...}}`. Inexact inputs yield only as many terms as the precision determines; the `Pi` example exposes the famous large term `292`. `ContinuedFraction` is Listable.
