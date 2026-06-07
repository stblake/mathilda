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

**Algorithm.** `Plus` carries `FLAT | ORDERLESS | LISTABLE | NUMERICFUNCTION | ONEIDENTITY`, so by the time `builtin_plus` runs the evaluator has already flattened nested `Plus`, threaded over `List` args, and canonically sorted the arguments. The builtin handles the numeric folding and like-term collection. `Plus[]` is 0, `Plus[x]` is x.

It first distributes `Times[-1, Plus[...]]` over the outer sum (`is_neg_of_plus`) so cancellations like `a + b - (a+b) -> 0` actually collapse. Inexact contagion (`numeric_contagion_args`) numericalises exact numeric parts in-place when any summand is an inexact Real/MPFR. The main pass splits each term into a numeric coefficient and a base via `get_coeff_base` (`c·base`, with bare numbers having a null base), then groups by structurally equal base (`TermGroup` array), summing coefficients within each group. Pure numeric terms fold into a single running total. The result is reassembled as `Plus` of `coeff·base` terms, dropping zero-coefficient groups.

**Data structures.** A linear `TermGroup[]` array of `(base, coeff)` pairs; coefficient arithmetic goes through the exact integer/rational/bigint paths and `eval_and_free` for symbolic recombination. Like-term lookup is linear scan with `expr_eq`.

**Complexity / limits.** Roughly `O(n^2)` worst case from the linear like-term grouping over n summands (after the evaluator's canonical sort, equal bases are typically already adjacent).

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
- Source: [`src/plus.c`](https://github.com/stblake/mathilda/blob/main/src/plus.c)
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
