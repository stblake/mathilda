# Rational

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Rational[n, d]`**

represents the rational number n/d.

<details>
<summary>Notes</summary>

When n and d are integers, Rational auto-reduces by gcd, normalises the sign onto the numerator, and collapses to an Integer when d == 1. Rationals propagate through Plus / Times exactly via GMP.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Rational[15, 5]
Out[1]= 3
```

### Applications (5)

```mathematica
In[2]:= Rational[6, 4]
Out[2]= 3/2

In[3]:= Rational[10, 2]
Out[3]= 5

In[4]:= Rational[-3, -9]
Out[4]= 1/3

In[5]:= 1/2 + 1/3 + 1/6
Out[5]= 1

In[6]:= Sum[1/k^2, {k, 1, 10}]
Out[6]= 1968329/1270080
```

## Implementation notes

`Rational[n, d]` is the internal head for exact rationals. `builtin_rational` only fires for two integer arguments: it calls `make_rational(n, d)` to reduce to lowest terms with a positive denominator. If the input is already in canonical form (no reduction happened) it returns `NULL` so the structural `Rational[n, d]` is left as-is; otherwise it returns the reduced form (an `EXPR_INTEGER` when the denominator becomes 1). Division by zero emits `Power::infy` and returns `ComplexInfinity` (or `Indeterminate` for `0/0`). Non-integer arguments return `NULL`.

- Automatically simplifies to lowest terms (e.g. `Rational[15, 5]` evaluates to `3`, `Rational[2, 4]` evaluates to `Rational[1, 2]`).
- Returns `Indeterminate` when `n` and `d` are both `0` (e.g. `Rational[0, 0]`).
- Returns `ComplexInfinity` when `n` is non-zero and `d` is `0` (e.g. `Rational[1, 0]`).

**Attributes:** `Protected`.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_arc_exact.c`](https://github.com/stblake/mathilda/blob/main/tests/test_arc_exact.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)

## Notes & additional examples

### Notes

`Rational[n, d]` represents the rational number `n/d`. With integer arguments it
reduces to lowest terms, moves the sign to the numerator, and becomes an
`Integer` when `d` divides `n`. The head of any non-integer fraction such as
`1/2` is `Rational`.
