# Equal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

lhs == rhs or Equal\[lhs, rhs\] tests mathematical equality. Numeric arguments decide directly (Integer / Rational exact comparison; Real / MPFR comparison with precision tolerance); structurally identical symbolic forms decide True; otherwise the call stays unevaluated as a symbolic equation. Equal threads over Lists pairwise; chained Equal becomes Inequality. Following IEEE 754 / ISO 60559, Indeterminate is unordered with every value including itself, so any Indeterminate argument gives False.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 2 == 2
Out[1]= True

In[2]:= 1 == 1.
Out[2]= True

In[3]:= a == b
Out[3]= a == b
```

## Implementation notes

**Algorithm.** `builtin_equal` walks adjacent argument pairs. For each pair it first tests structural identity (`expr_eq`); if that fails it calls `compare_numeric`. `compare_numeric` does exact GMP comparison (`mpz_cmp`) when both sides are integer-like (so `10^30 == 10^30 + 1` is correctly False even past 2^53), exact `long double` cross-multiplied comparison when both are rational/integer, and otherwise a tolerance comparison on the doubles (relative tolerance `2^-46`) so machine reals that agree to ~14 digits compare equal. A pair compares equal → continue; a decidable non-equal pair (or two distinct "raw data" leaves, via `is_raw_data`) → return `False` immediately. If some pair is undecidable (symbolic), the whole call returns NULL (unevaluated). All-equal returns `True`. `Equal[]`/`Equal[x]` return `True`.

**Data structures.** Operates directly on the `Expr` argument array; numeric extraction goes through `get_numeric_value` (double + exact rational num/den + exactness flag) and GMP `mpz_t` for big integers.

- Numeric arguments are compared by value, so `2 == 2.0` is `True`.
- For symbolic arguments that cannot be decided, the expression is returned
  unevaluated (`x == y`).
- `Equal` is `Orderless` for the equality test but preserves Mathematica's
  printed form.
- An `Indeterminate` argument gives `False`, per IEEE 754 — see
  [Indeterminate and IEEE unordered comparison](#indeterminate-and-ieee-unordered-comparison).

**Attributes:** `Protected`.

## See also

[Orderless](../../expression-information/Orderless/)

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_comparisons.c`](https://github.com/stblake/mathilda/blob/main/tests/test_comparisons.c)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_eliminate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eliminate.c)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)

## Notes & additional examples

### Notes

Unlike `SameQ`, `Equal` (`==`) tests mathematical equality, so `1 == 1.` is `True`. When equality cannot be decided, the call stays unevaluated as a symbolic equation.
