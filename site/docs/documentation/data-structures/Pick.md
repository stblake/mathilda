# Pick

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pick[expr, sel]`**

Picks out the elements of expr for which the corresponding element of sel is True.

**`Pick[expr, sel, patt]`**

Picks out the elements of expr for which the corresponding element of sel matches patt. Operates at all levels; sel must mirror the structure of expr, and the head of expr is preserved. Returns unevaluated if the structures disagree.

## Examples (6)

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

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
