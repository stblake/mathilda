# Equivalent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Equivalent[e1, e2, ...]`**

The logical equivalence e1 \\[Equivalent\] e2 \\[Equivalent\] ...: True when all of the ei have the same truth value.  Folds literal Booleans and cancels duplicate arguments; Equivalent\[\] and Equivalent\[e\] are True.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Equivalent[True, a, b]
Out[1]= a && b

In[2]:= Equivalent[False, a]
Out[2]= Not[a]

In[3]:= Equivalent[True, False]
Out[3]= False
```

### Applications (5)

All arguments share one truth value

```mathematica
In[4]:= Equivalent[True, True]
Out[4]= True
```

A True and a False together

```mathematica
In[5]:= Equivalent[True, False]
Out[5]= False
```

A literal True forces the remaining atoms

```mathematica
In[6]:= Equivalent[True, a, b]
Out[6]= a && b
```

A literal False negates the remaining atom

```mathematica
In[7]:= Equivalent[False, a]
Out[7]= Not[a]
```

Distinct symbolic atoms stay symbolic

```mathematica
In[8]:= Equivalent[p, q]
Out[8]= Equivalent[p, q]
```

## Implementation notes

**Algorithm.** `Equivalent` is `Flat, Orderless, OneIdentity, Protected`, so
nesting is flattened and arguments are sorted before `builtin_equivalent` runs.
`Equivalent[]` and `Equivalent[e]` are trivially `True`. Otherwise one pass
records whether a literal `True` and/or `False` is present and collects the
distinct non-literal atoms (duplicates dropped by `expr_eq`). A `True` together
with a `False` forces `False`. If a single truth value is present it forces every
remaining atom — `True` leaves each atom as-is, `False` wraps each in `Not` — and
the results are joined with `And` (a lone atom is returned directly). With no
literals but duplicates removed, a smaller `Equivalent` is rebuilt; distinct
symbolic atoms with nothing to simplify return `NULL` (stay symbolic).

**Data structures.** A single `malloc`'d array of borrowed pointers to the
distinct atoms, plus two `bool`s for the literal flags; only the surviving atoms
are deep-copied into the `And`/`Equivalent`/`Not` result.

**Complexity / limits.** The dedup check is pairwise (`O(n²)` `expr_eq`), adequate
for small boolean expressions. Simplification is structural; the full
truth-functional meaning is left to `LogicalExpand`/`Reduce`/`FindInstance`, which
expand `Equivalent[a1, …, an]` to the cyclic conjunction
`Implies[a1, a2] && … && Implies[an, a1]`.

**Attributes:** `Flat`, `OneIdentity`, `Orderless`, `Protected`.

## References

**See also:** [Flat](../../expression-information/Flat/), [Orderless](../../expression-information/Orderless/), [OneIdentity](../../expression-information/OneIdentity/), [LogicalExpand](../../solutions-of-equations/LogicalExpand/), [Reduce](../../solutions-of-equations/Reduce/), [FindInstance](../../solutions-of-equations/FindInstance/)

- Source: [`src/boolean.c`](https://github.com/stblake/mathilda/blob/main/src/boolean.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_boolean.c`](https://github.com/stblake/mathilda/blob/main/tests/test_boolean.c)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)

## Notes & additional examples

### Notes

`Equivalent[e1, e2, …]` is `True` when all of the `ei` share one truth value — all
`True` or all `False`. It is `Flat`, `Orderless` and `OneIdentity`; `Equivalent[]`
and `Equivalent[e]` are `True`.

Evaluation folds the literals and cancels duplicates: a `True` together with a
`False` gives `False`, a single literal forces each remaining atom (`True` leaves
it, `False` negates it) and the results are joined with `And`. `LogicalExpand`,
`Reduce` and `FindInstance` expand `Equivalent` to the cyclic conjunction
`Implies[a1, a2] && … && Implies[an, a1]`.
