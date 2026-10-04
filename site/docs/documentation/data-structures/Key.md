# Key

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Key[k]`**

Represents the key k of an association. Key\[k\]\[assoc\] gives the value at k (Missing\["KeyAbsent", k\] if absent); assoc\[\[Key\[k\]\]\], Lookup and the key-spec arguments of GroupBy, SortBy and JoinAcross accept it, and Key\[k\] names a key literally even when k is an integer or a list.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Key["b"][<|"a" -> 1, "b" -> 2|>]
Out[1]= 2

In[2]:= Lookup[<|1 -> "one", {1, 2} -> "pair"|>, Key[{1, 2}]]
Out[2]= "pair"

In[3]:= Key["z"][<|"a" -> 1|>]
Out[3]= Missing["KeyAbsent", "z"]
```

### Applications (3)

The curried lookup operator

```mathematica
In[4]:= Key["a"][<|"a" -> 1, "b" -> 2|>]
Out[4]= 1
```

Key as a Part specification

```mathematica
In[5]:= <|a -> 10, b -> 20|>[[Key[a]]]
Out[5]= 10
```

Name a field in a record pipeline

```mathematica
In[6]:= SortBy[{<|"n" -> 3|>, <|"n" -> 1|>, <|"n" -> 2|>}, Key["n"]]
Out[6]= {<|"n" -> 1|>, <|"n" -> 2|>, <|"n" -> 3|>}
```

## Implementation notes

**Algorithm.** `Key` is an inert wrapper: `Key[k]` represents the key `k` of an
association and does not evaluate on its own (it is registered `PROTECTED` with
no builtin body in `assoc_ops.c`). It acquires meaning in three places. The
operator form `Key[k][assoc]` is handled in the evaluator (`eval.c`): it looks
`k` up and returns the value or `Missing["KeyAbsent", k]` — the curried
complement of `assoc[Key[k]]`. As a `Part` specification, `assoc[[Key[k]]]`
reaches `k` through the association (`part.c`). And wherever a key spec is read —
`Lookup`, `KeyDrop`, `JoinAcross` — `Key[k]` is unwrapped to `k`.

**Data structures.** `Key[k]` is a one-argument `EXPR_FUNCTION`; the lookups it
drives go through `assoc_lookup_value`, i.e. the association's cached
open-addressing key index.

**Complexity / limits.** The wrapper itself is `O(1)`; the lookup it triggers is
`O(1)` amortised via the key index. `Key` exists so that record pipelines such as
`GroupBy[records, Key["field"]]` and `SortBy[records, Key["field"]]` can name a
key as an ordinary operator; it is distinct from the bare key only where a key is
itself something that would otherwise evaluate.

**Attributes:** `Protected`.

## References

**See also:** [GroupBy](../../data-structures/GroupBy/), [SortBy](../../data-structures/SortBy/), [JoinAcross](../../data-structures/JoinAcross/), [MapAt](../../data-structures/MapAt/), [Extract](../../structural-manipulation/Extract/)

- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)

## Notes & additional examples

### Notes

`Key[k]` names the key `k` of an association and is inert on its own. It acquires
meaning in three places: the operator form `Key[k][assoc]` returns the value at
`k` (or `Missing["KeyAbsent", k]`); as a part specification `assoc[[Key[k]]]` it
reaches that entry; and anywhere a key is read — `Lookup`, `KeyDrop`,
`JoinAcross` — `Key[k]` unwraps to `k`. Its purpose is to let a key act as an
ordinary operator, so record pipelines such as `SortBy[rows, Key["field"]]` and
`GroupBy[rows, Key["field"]]` can select a field by name. Wrapping is needed only
when the bare key would otherwise evaluate to something else.
