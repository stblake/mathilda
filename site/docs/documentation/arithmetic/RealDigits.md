# RealDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RealDigits[x] gives a list {digits, exp} of the digits in the approximate real number x together with the exponent such that the first digit is the coefficient of 10^(exp - 1).`**

**`RealDigits[x, b] gives base-b digits.`**

**`RealDigits[x, b, len] gives len digits.`**

**`RealDigits[x, b, len, n] gives len digits starting from the coefficient of b^n.`**

<details>
<summary>Notes</summary>

For rationals with non-terminating expansions the digit list ends in a nested list of the recurring block.  For inexact (machine or MPFR) reals, digits beyond the available precision are returned as Indeterminate.  The sign of x is discarded.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= RealDigits[123.55555]
Out[1]= {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}

In[2]:= RealDigits[Pi, 10, 25]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, 1}

In[3]:= RealDigits[19/7]
Out[3]= {{2, {7, 1, 4, 2, 8, 5}}, 1}

In[4]:= RealDigits[5.635, 10, 20]
Out[4]= {{5, 6, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}

In[5]:= RealDigits[Pi, 10, 20, -5]
Out[5]= {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3}, -4}

In[6]:= RealDigits[1.234, 2, 15]
Out[6]= {{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}
```

### Applications (6)

```mathematica
In[1]:= RealDigits[123.456]
Out[1]= {{1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 3}
```

For an exact rational with a non-terminating decimal expansion, the digit list
ends in a *nested* list giving the recurring block — here the period-6 cycle of
`1/7`, and the mixed `22/7 = 3.142857142857…`:

```mathematica
In[1]:= RealDigits[1/7]
Out[1]= {{{1, 4, 2, 8, 5, 7}}, 0}

In[2]:= RealDigits[22/7]
Out[2]= {{3, {1, 4, 2, 8, 5, 7}}, 1}
```

High-precision constants expose their digits directly. Thirty digits of `π`
from a 30-digit MPFR value, and the first ten significant base-10 digits of a
40-digit `π`:

```mathematica
In[1]:= RealDigits[N[Pi, 30]]
Out[1]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3, 8, 3, 2, 8}, 1}

In[2]:= RealDigits[N[Pi, 40], 10, 10]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3}, 1}
```

`RealDigits` also works in other bases — the binary expansion of `255` is eight
ones:

```mathematica
In[1]:= RealDigits[255, 2]
Out[1]= {{1, 1, 1, 1, 1, 1, 1, 1}, 8}
```

## Algorithm

real.c

RealDigits builtin -- positional-notation digit expansion.

```text
  RealDigits[x]               default base 10, length set by Precision[x].
  RealDigits[x, b]            base b, length set by Precision[x] / Log10[b].
  RealDigits[x, b, len]       exactly `len` digits, MSD-first.
  RealDigits[x, b, len, n]    `len` digits, first one = coefficient of b^n.

Result form is `{ digits-list, exp }`.  The first element of digits-list is
the coefficient of b^(exp-1).  Sign of x is discarded.  Exact rationals
```

with non-terminating base-b expansions return a list ending in a nested

```text
list of the recurring block.  Inexact (machine or MPFR) reals get
```

Indeterminate for any requested digit beyond the available precision.

```text
  x can be: Integer, BigInt, Rational[n,d], Real (machine), or
            EXPR_MPFR (arbitrary precision; USE_MPFR builds only).
```

The general algorithm scales |x| by base^(-low) where `low` is the lowest digit position we need, floors to an integer N, and reads off the base-b

```text
digits of N (padding with leading zeros as needed).  This single GMP /
MPFR shift handles every numeric type uniformly.  For the special case of
```

an exact rational with no explicit `len`, a remainder-tracked long division detects terminating vs recurring expansions and emits the nested-list form.

```text
Implementation only supports integer bases b >= 2.  Non-integer bases
```

(e.g. GoldenRatio) emit a `::ibase` diagnostic and leave the call unevaluated -- adding them requires MPFR floor-iteration and has been deferred.

## Implementation notes

**Algorithm.** `builtin_realdigits` returns the digit list of a real number, in the standard `{digits, exponent}` form. It accepts `RealDigits[x]`, `RealDigits[x, b]`, `RealDigits[x, b, len]`, `RealDigits[x, b, len, p]` (1–4 args; wrong count emits `RealDigits::argb`). `x` is classified by `rd_classify`; concrete non-real `Complex` input emits `RealDigits::realx`, and symbolic constants (Pi, E, …) are numericalised only once enough precision context (base and length) is known. The base defaults to 10, must be an integer ≥ 2 (`RealDigits::ibase` otherwise) and fit in `unsigned long`. Digits are extracted by repeated scaled-floor / MPFR shifting in the requested base, honouring the optional length and starting-position arguments.

**Data structures.** GMP `mpz_t` for the base and integer parts; MPFR for the fractional digit extraction when built. Output is a `List` of digits paired with an integer exponent.

- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under
  `USE_MPFR`) arbitrary-precision `MPFR` numbers. The sign of `x` is
  discarded.
- For integers and rationals with terminating base-`b` expansions, the
  digit list is flat. For rationals with non-terminating expansions, the
  list ends in a nested list giving the recurring block:
  `RealDigits[19/7]` returns `{{2, {7, 1, 4, 2, 8, 5}}, 1}`.
- For inexact (machine or MPFR) reals, the default `len` is
  `Round[Precision[x] / Log10[b]]`. Requesting more digits than the
  precision allows produces `Indeterminate` at the LSB end. The digits
  themselves use the canonical round-to-nearest representation supplied
  by MPFR -- so `RealDigits[123.55555]` returns the literal decimal
  digits, not the binary IEEE tail.
- Symbolic numeric constants such as `Pi`, `E`, and `GoldenRatio` are
  numericalized to MPFR at the matching precision when an explicit `len`
  is given. `RealDigits[Pi]` (no `len`) is left unevaluated.
- `RealDigits[0]` returns `{{0}, 0}`. `RealDigits[0.]` returns
  `{{0}, -Floor[Accuracy[0.]]}` — `{{0}, -323}` for machine precision
  (`Accuracy[0.] ≈ 323.607`), and `{{0}, -p}` for an MPFR zero
  `0``p` of precision `p` digits.
- Bases must be integers `>= 2`. Non-integer (Real / Rational) bases
  trigger `RealDigits::ibase` and leave the call unevaluated.
  Non-integer-base expansions (e.g. `GoldenRatio`) are not yet supported.
- `FromDigits` can be used as the inverse of `RealDigits` for the
  integer / terminating-rational case.

**Attributes:** `Listable`, `Protected`.

## See also

[Rational](../../arithmetic/Rational/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [GoldenRatio](../../mathematical-constants/GoldenRatio/), [FromDigits](../../arithmetic/FromDigits/)

## References

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_real_digits.c`](https://github.com/stblake/mathilda/blob/main/tests/test_real_digits.c)

## Notes & additional examples

### Notes

`RealDigits[x]` gives `{digits, exp}` where the first digit is the coefficient of
`10^(exp - 1)`. `RealDigits[x, b]` uses base `b`; `RealDigits[x, b, len]` returns
`len` digits; `RealDigits[x, b, len, n]` starts from the coefficient of `b^n`.
For rationals with non-terminating expansions the digit list ends in a nested
list of the recurring block. For inexact reals, digits beyond the available
precision are returned as `Indeterminate`. The sign of `x` is discarded.
