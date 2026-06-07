# Plus

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x + y + ... or Plus[x, y, ...] represents a sum of terms.
Plus is Flat, Orderless, OneIdentity, Listable, and NumericFunction:
nested Plus is auto-flattened, terms are sorted into canonical order,
like terms are combined, and integer arguments are summed exactly via
GMP (with int64 fast path and BigInt overflow promotion).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3 x
```

## Implementation notes

- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and collects like terms (e.g., `x + 2x` becomes `3x`).
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on arbitrary-precision integer addition.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on normal forms for sums.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2^100 + 3^50
Out[1]= 1267651318126217093349291975625
```

```mathematica
In[1]:= 1/2 + 1/3 + 1/6
Out[1]= 1
```

```mathematica
In[1]:= a + a + b + 2 a
Out[1]= 4 a + b
```

### Notes

Plus is `Flat`, `Orderless`, and `OneIdentity`, so nested sums are flattened and
terms are sorted into canonical order before combination. Integer addition
auto-promotes to GMP bigints on overflow, as in `2^100 + 3^50`. Exact rationals
are added with a common denominator and the result is reduced to lowest terms,
collapsing to an integer when the denominator divides out. Like terms are
collected by their symbolic factor, so `a + a + 2 a` becomes `4 a`.
