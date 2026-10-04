# Pick

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pick[expr, sel]`**

Picks out the elements of expr for which the corresponding element of sel is True.

**`Pick[expr, sel, patt]`**

Picks out the elements of expr for which the corresponding element of sel matches patt. Operates at all levels; sel must mirror the structure of expr, and the head of expr is preserved. Returns unevaluated if the structures disagree.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Catenate[{<|"a" -> 1|>, <|"b" -> 2|>}]
Out[1]= {1, 2}

In[2]:= Insert[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "d" -> 4, Key["b"]]
Out[2]= <|"a" -> 1, "d" -> 4, "b" -> 2, "c" -> 3|>

In[3]:= Insert[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "a" -> 9, 3]
Out[3]= <|"b" -> 2, "a" -> 9, "c" -> 3|>

In[4]:= Pick[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {True, False, True}]
Out[4]= <|"a" -> 1, "c" -> 3|>

In[5]:= Partition[<|"a" -> 1, "b" -> 2|>, 1]
Out[5]= Partition[<|"a" -> 1, "b" -> 2|>, 1]

In[6]:= Reverse[<|"x" -> {1, 2}, "y" -> {3, 4}|>, 2]
Out[6]= <|"x" -> {2, 1}, "y" -> {4, 3}|>
```

### Applications (3)

Two-argument form: keep where True

```mathematica
In[7]:= Pick[{a, b, c, d}, {True, False, True, False}]
Out[7]= {a, c}
```

Keep where the selector matches 1

```mathematica
In[8]:= Pick[{1, 2, 3, 4, 5}, {1, 0, 1, 0, 1}, 1]
Out[8]= {1, 3, 5}
```

Selector matched against a pattern

```mathematica
In[9]:= Pick[{1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}, _?EvenQ]
Out[9]= {2, 4, 6}
```

## Algorithm

Pick[expr, sel] / Pick[expr, sel, patt] — select elements of `expr` whose positionally-corresponding element of `sel` matches `patt` (literal `True` for the two-argument form).

The selector array mirrors the structure of `expr`, so the walk is a simultaneous recursive descent over both trees. At each element:

```text
  - selector matches `patt`            -> keep the whole `expr` element,
  - selector is compound and no match  -> recurse into the pair,
  - selector is atomic and no match    -> drop the element.
```

The head at every level comes from `expr`, never from `sel`, so Pick[f[a, b, c], {True, False, True}] is f[a, c].

Any structural disagreement between the two trees (differing lengths, or a compound selector against an atomic expression element) is not an error: the builtin returns NULL and the evaluator leaves `Pick[...]` unevaluated, matching Mathematica's Pick::incomp behaviour of returning the input. Detection is exact — a mismatch found arbitrarily deep aborts the whole call rather than yielding a partially picked result.

## Implementation notes

**Algorithm.** `Pick[expr, sel]` / `Pick[expr, sel, patt]` keeps the elements of
`expr` whose positionally-corresponding element of the selector `sel` matches `patt`
(literal `True` for the two-argument form). The selector mirrors `expr`'s structure,
so `pick_rec` descends both trees together: at each element, a selector that matches
`patt` keeps the whole `expr` element; a *compound* selector that did not match as a
whole recurses into the pair; an *atomic* non-matching selector drops the element.
The head at every level comes from `expr`, never `sel`. An association picks by
position over its entries — the value is compared and the entry (key included) kept;
an association used *as* a selector is atomic.

**Data structures.** The matcher (`match` via `pick_selects`) with a throwaway
`MatchEnv`; a per-level `kept` array rebuilt into a node with `expr`'s head.

**Complexity / limits.** O(size) with a match test per element. Any structural
disagreement — differing lengths, or a compound selector against an atomic
expression — is detected exactly (even arbitrarily deep) and aborts the whole call to
`NULL`, so `Pick[…]` is returned unevaluated, matching Mathematica's `Pick::incomp`.

- `Insert` removes an existing entry with the inserted key, so the new position
  wins; a non-rule element returns the association unchanged, and an
  out-of-range position or absent key leaves the call unevaluated
  (Mathematica 15).
- An association used as a `Pick` selector is atomic: the whole expression if
  it matches the pattern, else `Sequence[]`, as in Mathematica 15 (it no longer
  builds malformed `Rule[]` nodes).

**Attributes:** `Protected`.

## References

**See also:** [Catenate](../../data-structures/Catenate/), [Insert](../../data-structures/Insert/), [Partition](../../data-structures/Partition/)

- Source: [`src/list/pick.c`](https://github.com/stblake/mathilda/blob/main/src/list/pick.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)

## Notes & additional examples

### Notes

The selector must mirror the structure of the first argument; `Pick` walks both in
step and keeps each element whose selector matches the pattern (`True` by default).
The result keeps the first argument's head. A shape mismatch leaves `Pick[…]`
unevaluated rather than erroring.
