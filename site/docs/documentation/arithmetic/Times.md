# Times

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x \* y \* ... or Times\[x, y, ...\] represents a product of terms. Times is Flat, Orderless, OneIdentity, Listable, and NumericFunction: nested Times is auto-flattened, factors are sorted, like factors are merged into Power, and integer products use exact GMP arithmetic. Numeric zero collapses the product; a Plus factor is left distributed (use Expand to distribute).

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= 2 * x * 3 * y
Out[1]= 6 x y

In[2]:= I * I
Out[2]= -1

In[3]:= 14/Sqrt[10]
Out[3]= 7 Sqrt[2/5]
```

### Applications (3)

```mathematica
In[4]:= 2^64 * 3
Out[4]= 55340232221128654848

In[5]:= (1/3)*(3/7)*7
Out[5]= 1

In[6]:= x y x z
Out[6]= x^2 y z
```

## Implementation notes

**Algorithm.** `Times` carries `FLAT | ORDERLESS | LISTABLE | NUMERICFUNCTION | ONEIDENTITY`, so the evaluator flattens nested `Times`, threads over `List` args, and sorts canonically before `builtin_times` runs. The builtin does numeric folding and exponent collection.

It accumulates the rational/numeric coefficient into a running `num_prod` (with inexact contagion when a Real/MPFR factor is present); a numeric factor that cannot fold (e.g. an integer carrying a pending radical) is stashed as its own group. Each remaining factor is split into a `(base, exponent)` pair (a `BasePower`; bare factors get exponent 1, `Power[b,e]` is unpacked), and factors with structurally equal base are merged by summing exponents (`Plus` over the two exponents). A radical-canonicalisation pass folds factors of a positive integer base `b ≥ 2` appearing in the accumulated rational coefficient into a `Power[b, q]` group with rational q. The result is rebuilt as the numeric coefficient times the surviving `Power` groups; `0` short-circuits and a zero result returns 0.

**Data structures.** The `BasePower` struct (`{base, exponent}` — shared with `trig_canon.h`) in a linear `groups[]` array; coefficient arithmetic on exact int/rational/bigint with `eval_and_free` for symbolic recombination. Like-base lookup is linear scan with `expr_eq`. An overflow guard returns an `Overflow[]` sentinel.

**Complexity / limits.** Roughly `O(n^2)` from the linear base-grouping over n factors (canonical sort tends to make equal bases adjacent).

- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.
- **Sqrt-coefficient absorption**: a non-trivial rational/integer coefficient
  combined with a single `Power[r, -1/2]` group folds into the canonical
  Mathematica form `sign(c) * Sqrt[c^2 / r]`, then Power's existing
  rational-base extraction pulls out any newly-exposed perfect square:
  `14/Sqrt[10] → 7 Sqrt[2/5]`, `78/Sqrt[66] → 13 Sqrt[6/11]`,
  `(3/5)/Sqrt[2/5] → 3/Sqrt[10]`. Skipped when no square is actually
  extracted (e.g. `2/Sqrt[30]` stays as-is to preserve like-term
  collection in `Plus`).

**Attributes:** `Flat`, `Listable`, `NumericFunction`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [Power](../../arithmetic/Power/), [I](../../mathematical-constants/I/), [Plus](../../arithmetic/Plus/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on arbitrary-precision integer multiplication.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on canonical forms for products.
- Source: [`src/times.c`](https://github.com/stblake/mathilda/blob/main/src/times.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_arc_exact.c`](https://github.com/stblake/mathilda/blob/main/tests/test_arc_exact.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_chop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chop.c)

## Notes & additional examples

### Notes

Times shares the `Flat`, `Orderless`, and `OneIdentity` attributes with Plus, so
factors are flattened and canonically ordered. Repeated symbolic factors are
gathered into powers, turning `x y x z` into `x^2 y z`. Rational factors are
multiplied and reduced to lowest terms, and a product of reciprocals can cancel
exactly to `1`. Integer products promote to GMP bigints automatically when they
exceed machine-word range.
