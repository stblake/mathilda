# Rational

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Rational[n, d]
    represents the rational number n/d.
When n and d are integers, Rational auto-reduces by gcd, normalises
the sign onto the numerator, and collapses to an Integer when d == 1.
Rationals propagate through Plus / Times exactly via GMP.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Rational[15, 5]
Out[1]= 3
```

## Implementation notes

`Rational[n, d]` is the internal head for exact rationals. `builtin_rational` only fires for two integer arguments: it calls `make_rational(n, d)` to reduce to lowest terms with a positive denominator. If the input is already in canonical form (no reduction happened) it returns `NULL` so the structural `Rational[n, d]` is left as-is; otherwise it returns the reduced form (an `EXPR_INTEGER` when the denominator becomes 1). Division by zero emits `Power::infy` and returns `ComplexInfinity` (or `Indeterminate` for `0/0`). Non-integer arguments return `NULL`.

- Automatically simplifies to lowest terms (e.g. `Rational[15, 5]` evaluates to `3`, `Rational[2, 4]` evaluates to `Rational[1, 2]`).
- Returns `Indeterminate` when `n` and `d` are both `0` (e.g. `Rational[0, 0]`).
- Returns `ComplexInfinity` when `n` is non-zero and `d` is `0` (e.g. `Rational[1, 0]`).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Rational[6, 4]
Out[1]= 3/2
```

`Rational` auto-reduces by the gcd, normalises the sign onto the numerator, and
collapses to an `Integer` whenever the denominator divides the numerator:

```mathematica
In[1]:= Rational[10, 2]
Out[1]= 5

In[2]:= Rational[-3, -9]
Out[2]= 1/3
```

Because rationals propagate exactly through `Plus` and `Times` via GMP, exact
sums never drift into floating point — a partial sum of the Basel series stays a
single reduced fraction:

```mathematica
In[1]:= 1/2 + 1/3 + 1/6
Out[1]= 1

In[2]:= Sum[1/k^2, {k, 1, 10}]
Out[2]= 1968329/1270080
```

### Notes

`Rational[n, d]` represents the rational number `n/d`. With integer arguments it
reduces to lowest terms, moves the sign to the numerator, and becomes an
`Integer` when `d` divides `n`. The head of any non-integer fraction such as
`1/2` is `Rational`.
