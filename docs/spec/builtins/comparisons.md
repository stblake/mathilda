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

Note this differs from the reference, where `Indeterminate == Indeterminate` is
`True`.

## Equal
Tests whether two (or more) expressions are mathematically equal.
- `Equal[x, y]` (written `x == y`): `True` if `x` and `y` are equal, `False` if
  they are provably unequal, otherwise stays symbolic.
- `Equal[x, y, z, ...]` (written `x == y == z`): tests that all arguments are
  equal.

**Features**:
- Numeric arguments are compared by value, so `2 == 2.0` is `True`.
- Approximate numbers are equal when they agree to within the uncertainty of
  the **less accurate** one, ignoring its last 7 bits (about 2.1 decimal
  digits). With precision `p` bits (53 for a machine real; an approximate zero
  has none) and accuracy `a = p - log2|x|` bits (a zero: its own bit
  precision), the operand of least accuracy governs and `x == y` holds iff
  `|x - y| <= 2^(k - a)`, `k = 7 min(p, 53)/53`. An exact operand (`1`, `Pi`,
  `Sqrt[2]`) has unlimited accuracy, so it compares at the approximate side's
  precision. Two machine reals keep the relative `2^-46` rule. Complex numbers
  compare part by part.
  `N[1 - 10^-20, 30] == 1` is `False` and `N[1 - 10^-29, 30] == 1` is `True`.
  Before this, an arbitrary-precision operand was read as a double and compared
  at the machine tolerance, so `N[1 - 1192/10^30, 30] == 1` was `True` and
  `N[1 - 1192/10^30, 30] < 1` was `False`.
  `N[x, d]` has a precision of `ceil(d log2 10)` bits (`Precision[N[1, 30]]`
  is 30.; the mantissa's 16 guard bits do not count), so the window at 1 is
  `2^(7 - 100) ~ 1.0 10^-28` at 30 digits.
- The ordering heads agree with it: inside the window `x < y` is `False` and
  `x <= y` is `True`; outside it the sign of `x - y` decides. The window holds
  beyond the double range too (`N[10^3000, 30] == N[10^3000 (1 + 10^-29), 30]`).
- The tolerance is not transitive, so a chain `a == b == c` needs **every**
  pair equal, not just neighbours: with `a = N[1, 30]`, `b = N[1 - 90/10^30, 30]`,
  `c = N[1 - 180/10^30, 30]`, `a == b` and `b == c` are `True`, `a == b == c`
  is `False`. Machine reals follow the same rule
  (`1. == 1. + 10.^-14 == 1. + 2 10.^-14` is `False`).
- An exact cancellation leaves an approximate zero as accurate as its
  operands: `N[1/34, 30] - N[1/34, 30]` has accuracy above 31 digits, not 30,
  so it is `< 10^-31`. The zero keeps that accuracy through `N[z, p]`
  (an approximate zero has no precision to lower: `N[0``80, 20]` is `0``80`),
  and a product or quotient scales it: `0``80 * 10^50` is `0``30` (known to
  within `10^-30`, so not `< 10^-40`), `0``80 / 10^35` is `0``115`.
- Two exact rationals compare exactly, whatever the size of their parts
  (v0.477): `1 == 1 + 10^-20` is `False`. It used to reach the tolerant
  machine comparison and test `True`.
- Closed-form numeric constants that carry no real ordering are decided by the
  exact zero-test: `I == 0` is `False`, `(1 + I) == (1 - I)` is `False`, and a
  complex radical such as `(-1 + (1/8)(-1 - I Sqrt[3])^3)/(-1 + (1/2)(-1 - I Sqrt[3]))
  == 0` is `True`. The test only fires when **both** sides are numeric
  (`NumericQ`), so a free symbol still stays symbolic (`x == 0`).
- For symbolic arguments that cannot be decided, the expression is returned
  unevaluated (`x == y`).
- A side is numericalised for the value comparison only when it is a numeric
  quantity (`NumericQ`), as in the reference (v0.6129). `N` maps into the
  arguments of an inert call, so `ff[5] == x` used to evaluate `ff[5.]` and
  fire a definition on real arguments; FeynCalc's `GA[5.]` guard printed
  `GA::noint` inside every `... == x` test. `Unequal`, `Less` and the other
  comparisons share the same fallback.
- Two associations are equal when their entries agree in order, with identical
  keys and `Equal` values: `<|a -> 1|> == <|a -> 1.|>` is `True`,
  `<|a -> 1, b -> 2|> == <|b -> 2, a -> 1|>` is `False`. `Unequal` uses the
  same test (see
  [Association atomicity](data-structures.md#association-atomicity-structural-functions)).
- Two Lists are compared element by element with `Equal` itself, recursively
  through nested lists, so the approximate-number tolerance applies inside them:
  `{0.1 + 0.2, 1.} == {0.3, 1.}` is `True`. Lists of different lengths are
  `False`, one decidably unequal position makes the whole comparison `False`
  (`{x, 1} == {y, 2}`), and otherwise it stays symbolic (`{x} == {y}`).
  `Unequal` uses the same test.
- `Equal` is `Orderless` for the equality test but preserves the reference's
  printed form.

```mathematica
In[1]:= {0.1 + 0.2} == {0.3}
Out[1]= True

In[2]:= {1, 2} == {1, 2, 3}
Out[2]= False

In[3]:= {x, 1} == {y, 2}
Out[3]= False
```
- An `Indeterminate` argument gives `False`, per IEEE 754 — see
  [Indeterminate and IEEE unordered comparison](#indeterminate-and-ieee-unordered-comparison).

**v0.6034.** A visible `NDArray[...]` operand compares through its nested List, as a packed list does: `NDArray[{3., 2.}] == {3., 2.}` is `True`, and arrays that differ in an element are `False`. `Unequal` likewise.

## Unequal
Tests whether expressions are unequal.
- `Unequal[x, y]` (written `x != y`): `True` if `x` and `y` are provably
  unequal, `False` if equal, otherwise stays symbolic.
- `Unequal[x, y, z, ...]` (written `x != y != z`): `True` only if all
  arguments are pairwise distinct. A run of `!=` parses to this one variadic
  call, as in the reference (`1 != 2 != 1` is `False`); mixed with an ordering
  operator it joins an `Inequality` instead (`a < b != c`).

**Features**:
- A pair containing `Indeterminate` counts as unequal, so
  `Indeterminate != Indeterminate` is `True`. An equal non-`Indeterminate` pair
  still decides the whole call `False`.
- Like `Equal`, an undecided pair of numeric constants is settled by the exact
  zero-test, so `I != 0` and `2 I Pi != 0` are `True`, while a free symbol stays
  symbolic (`x != 0`).

## Less
Tests strict ascending order.
- `Less[x, y]` (written `x < y`): `True` if `x` is less than `y`.
- `Less[x, y, z, ...]` (written `x < y < z`): tests the whole chain.
- An `Indeterminate` argument gives `False`; likewise for `LessEqual`,
  `Greater` and `GreaterEqual`.
- `Infinity` and `-Infinity` (and `DirectedInfinity[d]` for a real `d`) order
  against every real number and each other, for all four orderings:
  `0 < Infinity`, `-Infinity < Pi` and `Infinity <= Infinity` are `True`, and
  `Infinity < Infinity` is `False`. A symbolic or complex other side, or
  `ComplexInfinity`, stays unevaluated.

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
- `SameQ[x, y, z, ...]` (written `x === y === z`): `True` if all arguments are
  identical. A run of `===` parses to this one variadic call, as in
  the reference, so `1 === 1 === 1` is `True`; `(a === b) === c` stays nested.

**Features**:
- Unlike `Equal`, `SameQ` never stays symbolic and does not coerce numeric
  types, so `2 === 2.0` is `False`.
- Also unlike `Equal`, `SameQ` compares structure rather than value, so
  `Indeterminate === Indeterminate` is `True`.

## UnsameQ
Tests structural non-identity.
- `UnsameQ[x, y]` (written `x =!= y`): the logical negation of `SameQ`.
- `UnsameQ[x, y, z, ...]` (written `x =!= y =!= z`): `True` if no two
  arguments are identical.

```mathematica
In[1]:= 1 === 1 === 1
Out[1]= True

In[2]:= FullForm[Hold[a =!= b =!= c]]
Out[2]= Hold[UnsameQ[a, b, c]]

In[3]:= 1 =!= 2 =!= 1
Out[3]= False
```

## Inequality
The canonical internal form for chained comparisons.
- `Inequality[e1, op1, e2, op2, e3, ...]`: produced by the parser for
  expressions such as `a < b <= c`. Each `opk` is one of `Less`, `LessEqual`,
  `Greater`, `GreaterEqual`, `Equal` or `Unequal` (the last only when mixed with
  the others, `a < b != c`). The chain is `True` only when every adjacent
  comparison holds; an `Equal`/`Unequal` pair is decided by that builtin, list
  tolerance and exact zero-test included.

**Features**:
- Adjacent comparisons that evaluate to `True` are dropped; the result collapses
  to `True`, `False`, or a reduced `Inequality` over the undecided operands.
- An `Indeterminate` operand makes its adjacent pairs `False`, so the whole
  chain is `False`.
