# Comparisons

Relational operators compare expressions and return `True`, `False`, or stay
unevaluated (symbolic) when the relation cannot be decided. Equality tests come
in two flavours: mathematical equality (`Equal`, `Unequal`), which may stay
symbolic, and structural identity (`SameQ`, `UnsameQ`), which always returns a
boolean. Ordering operators (`Less`, `Greater`, and their `…Equal` variants)
support chained comparisons such as `a < b < c`, which the parser folds into an
`Inequality` expression.

## Indeterminate and IEEE unordered comparison

`Indeterminate` is Mathilda's NaN — it is what `0/0`, `Infinity - Infinity` and
`0 ComplexInfinity` normalise to. Following IEEE 754 / ISO 60559, it is
*unordered* with respect to every value, itself included:

| Comparison | Result with an `Indeterminate` operand |
|------------|----------------------------------------|
| `Equal`, `Less`, `LessEqual`, `Greater`, `GreaterEqual` | `False` |
| `Unequal` | `True` |
| `SameQ`, `UnsameQ` | unaffected — structural, see below |

```
Indeterminate == Indeterminate   (* False *)
Indeterminate != Indeterminate   (* True  *)
Indeterminate <  Indeterminate   (* False *)
Indeterminate == 1               (* False *)
0/0 == 0/0                       (* False *)
1 < Indeterminate < 3            (* False *)
```

Every argument of a chained comparison takes part in at least one adjacent
pair, so a single `Indeterminate` anywhere in the chain decides the whole
expression.

`SameQ` and `UnsameQ` are deliberately exempt: they test structural identity of
expressions rather than mathematical equality of values — the same distinction
IEEE draws between its comparison predicates and `totalOrder`, which does order
a NaN with itself. So `Indeterminate === Indeterminate` is `True`. This applies
to `Indeterminate` as an *operand*; a compound operand that merely contains one
(e.g. `{1, Indeterminate} == {1, Indeterminate}`) still decides by structure.

Note this differs from Mathematica, where `Indeterminate == Indeterminate` is
`True`.

## Equal
Tests whether two (or more) expressions are mathematically equal.
- `Equal[x, y]` (written `x == y`): `True` if `x` and `y` are equal, `False` if
  they are provably unequal, otherwise stays symbolic.
- `Equal[x, y, z, ...]` (written `x == y == z`): tests that all arguments are
  equal.

**Features**:
- Numeric arguments are compared by value, so `2 == 2.0` is `True`.
- For symbolic arguments that cannot be decided, the expression is returned
  unevaluated (`x == y`).
- `Equal` is `Orderless` for the equality test but preserves Mathematica's
  printed form.
- An `Indeterminate` argument gives `False`, per IEEE 754 — see
  [Indeterminate and IEEE unordered comparison](#indeterminate-and-ieee-unordered-comparison).

## Unequal
Tests whether expressions are unequal.
- `Unequal[x, y]` (written `x != y`): `True` if `x` and `y` are provably
  unequal, `False` if equal, otherwise stays symbolic.
- `Unequal[x, y, z, ...]`: `True` only if all arguments are pairwise distinct.

**Features**:
- A pair containing `Indeterminate` counts as unequal, so
  `Indeterminate != Indeterminate` is `True`. An equal non-`Indeterminate` pair
  still decides the whole call `False`.

## Less
Tests strict ascending order.
- `Less[x, y]` (written `x < y`): `True` if `x` is less than `y`.
- `Less[x, y, z, ...]` (written `x < y < z`): tests the whole chain.
- An `Indeterminate` argument gives `False`; likewise for `LessEqual`,
  `Greater` and `GreaterEqual`.

## LessEqual
Tests non-strict ascending order.
- `LessEqual[x, y]` (written `x <= y`): `True` if `x` is less than or equal to
  `y`.

## Greater
Tests strict descending order.
- `Greater[x, y]` (written `x > y`): `True` if `x` is greater than `y`.
- `Greater[x, y, z, ...]` (written `x > y > z`): tests the whole chain.

## GreaterEqual
Tests non-strict descending order.
- `GreaterEqual[x, y]` (written `x >= y`): `True` if `x` is greater than or
  equal to `y`.

## SameQ
Tests structural identity.
- `SameQ[x, y]` (written `x === y`): `True` if `x` and `y` are structurally
  identical, `False` otherwise. Always returns a boolean.

**Features**:
- Unlike `Equal`, `SameQ` never stays symbolic and does not coerce numeric
  types, so `2 === 2.0` is `False`.
- Also unlike `Equal`, `SameQ` compares structure rather than value, so
  `Indeterminate === Indeterminate` is `True`.

## UnsameQ
Tests structural non-identity.
- `UnsameQ[x, y]` (written `x =!= y`): the logical negation of `SameQ`.

## Inequality
The canonical internal form for chained comparisons.
- `Inequality[e1, op1, e2, op2, e3, ...]`: produced by the parser for
  expressions such as `a < b <= c`. Each `opk` is one of the relational heads
  above. The chain is `True` only when every adjacent comparison holds.

**Features**:
- Adjacent comparisons that evaluate to `True` are dropped; the result collapses
  to `True`, `False`, or a reduced `Inequality` over the undecided operands.
- An `Indeterminate` operand makes its adjacent pairs `False`, so the whole
  chain is `False`.
