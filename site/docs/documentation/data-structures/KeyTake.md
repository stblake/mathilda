# KeyTake

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyTake[assoc, {k1, ...}]`**

Gives the association of only the specified keys, in the requested order.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyTake[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"c", "a"}]
Out[1]= <|"c" -> 3, "a" -> 1|>

In[2]:= KeyTake[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, {"a"}]
Out[2]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

### Applications (1)

Keep only a and c

```mathematica
In[3]:= KeyTake[<|a -> 1, b -> 2, c -> 3|>, {a, c}]
Out[3]= <|a -> 1, c -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_keytake` is `key_drop_take(res, take=true)`. For a single
association it calls `assoc_key_select(assoc, karg, take=true)`, which keeps only the
listed keys while preserving the association's order; `KeyDrop` is the same routine
with `take=false`. The key argument may be one key or a `List` of keys. When the
first argument is a non-empty list of associations the call threads, delegating to
itself per element (the column-of-records form).

**Data structures.** The requested keys are loaded into a transient open-addressing
`KeyIndex` once, so each entry is classified by a single O(1) probe. The result is a
fresh `Association`; the compiled evaluator (B3) calls `assoc_key_select` directly,
with no call-node round-trip.

**Complexity / limits.** O(n + k) for `n` entries and `k` requested keys. A first
argument that is neither an association nor a threadable list leaves the call
unevaluated.

- Matches Mathematica 15: absent keys are skipped, and a key requested more
  than once is placed at its last occurrence (`{"c", "a", "c"}` gives
  `<|"a" -> .., "c" -> ..|>`).
- `O(n + m)`: the association's entries are hash-indexed once; `RuleDelayed`
  entries stay delayed.

**Attributes:** `Protected`.

## References

**See also:** [RuleDelayed](../../assignment-and-rules/RuleDelayed/)

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

The kept entries stay in the association's own order, not the order of the key list.
`KeyDrop` is the complement (remove the listed keys). Both thread over a list of
associations, which is the column-of-records form.
