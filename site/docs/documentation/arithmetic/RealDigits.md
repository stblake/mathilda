# RealDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RealDigits[x] gives a list {digits, exp} of the digits in the approximate real number x together with the exponent such that the first digit is the coefficient of 10^(exp - 1).
RealDigits[x, b] gives base-b digits.
RealDigits[x, b, len] gives len digits.
RealDigits[x, b, len, n] gives len digits starting from the coefficient of b^n.
For rationals with non-terminating expansions the digit list ends in a nested list of the recurring block.  For inexact (machine or MPFR) reals, digits beyond the available precision are returned as Indeterminate.  The sign of x is discarded.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RealDigits[123.55555]
Out[1]= {{1, 2, 3, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0}, 3}

In[2]:= RealDigits[Pi, 10, 25]
Out[2]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}

In[3]:= RealDigits[19/7]
Out[3]= {{2, {7, 1, 4, 2, 8, 5}}, 1}

In[4]:= RealDigits[5.635, 10, 20]
Out[4]= {{5, 6, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, 1}

In[5]:= RealDigits[Pi, 10, 20, -5]
Out[5]= {{9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate, Indeterminate}, -4}

In[6]:= RealDigits[1.234, 2, 15]
Out[6]= {{1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 1}
```

## Implementation notes

**Algorithm.** `builtin_realdigits` returns the digit list of a real number, in the Mathematica `{digits, exponent}` form. It accepts `RealDigits[x]`, `RealDigits[x, b]`, `RealDigits[x, b, len]`, `RealDigits[x, b, len, p]` (1–4 args; wrong count emits `RealDigits::argb`). `x` is classified by `rd_classify`; concrete non-real `Complex` input emits `RealDigits::realx`, and symbolic constants (Pi, E, …) are numericalised only once enough precision context (base and length) is known. The base defaults to 10, must be an integer ≥ 2 (`RealDigits::ibase` otherwise) and fit in `unsigned long`. Digits are extracted by repeated scaled-floor / MPFR shifting in the requested base, honouring the optional length and starting-position arguments.

**Data structures.** GMP `mpz_t` for the base and integer parts; MPFR for the fractional digit extraction when built. Output is a `List` of digits paired with an integer exponent.

- `Protected`, `Listable`. Threads over lists in any argument position.
- Works for `Integer`, `BigInt`, `Rational`, machine `Real`, and (under

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/real.c`](https://github.com/stblake/mathilda/blob/main/src/real.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
