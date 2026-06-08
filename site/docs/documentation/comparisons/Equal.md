# Equal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
lhs == rhs or Equal[lhs, rhs]
    tests mathematical equality. Numeric arguments decide directly
    (Integer / Rational exact comparison; Real / MPFR comparison with
    precision tolerance); structurally identical symbolic forms decide
    True; otherwise the call stays unevaluated as a symbolic equation.
Equal threads over Lists pairwise; chained Equal becomes Inequality.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_equal` walks adjacent argument pairs. For each pair it first tests structural identity (`expr_eq`); if that fails it calls `compare_numeric`. `compare_numeric` does exact GMP comparison (`mpz_cmp`) when both sides are integer-like (so `10^30 == 10^30 + 1` is correctly False even past 2^53), exact `long double` cross-multiplied comparison when both are rational/integer, and otherwise a tolerance comparison on the doubles (relative tolerance `2^-46`) so machine reals that agree to ~14 digits compare equal. A pair compares equal → continue; a decidable non-equal pair (or two distinct "raw data" leaves, via `is_raw_data`) → return `False` immediately. If some pair is undecidable (symbolic), the whole call returns NULL (unevaluated). All-equal returns `True`. `Equal[]`/`Equal[x]` return `True`.

**Data structures.** Operates directly on the `Expr` argument array; numeric extraction goes through `get_numeric_value` (double + exact rational num/den + exactness flag) and GMP `mpz_t` for big integers.

- Numeric arguments are compared by value, so `2 == 2.0` is `True`.
- For symbolic arguments that cannot be decided, the expression is returned

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 2 == 2
Out[1]= True

In[2]:= 1 == 1.
Out[2]= True

In[3]:= a == b
Out[3]= a == b
```

### Notes

Unlike `SameQ`, `Equal` (`==`) tests mathematical equality, so `1 == 1.` is `True`. When equality cannot be decided, the call stays unevaluated as a symbolic equation.
