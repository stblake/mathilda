# KeyIntersection

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyIntersection[{assoc1, assoc2, ...}]`**

Gives the list of associations restricted to the keys common to all of them, each in the key order of assoc1. Elements may also be rules or lists of rules.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>}]
Out[1]= {<|"b" -> 2|>, <|"b" -> 3|>}

In[2]:= KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "a" -> 4|>}]
Out[2]= {<|"a" -> 1, "b" -> 2|>, <|"a" -> 4, "b" -> 3|>}
```

### Applications (1)

```mathematica
In[3]:= KeyIntersection[{<|a -> 1, b -> 2, c -> 3|>, <|b -> 20, c -> 30, d -> 40|>}]
Out[3]= {<|b -> 2, c -> 3|>, <|b -> 20, c -> 30|>}
```

## Implementation notes

**Algorithm.** `builtin_keyintersection` restricts every association in a list to the
keys common to all of them. `as_assoc_array` first normalises the argument (its
elements may be associations, rules or lists of rules). It then scans the first
association's keys in order, keeping a key only when `assoc_lookup_value` finds it in
every other association; those common keys (in the first association's order) become
the key set each input is rebuilt against.

**Data structures.** Each association's persistent hash index backs the O(1)
membership probes; an `ExprBuf` holds the common keys, and one `ExprBuf` per
association accumulates the restricted entries.

**Complexity / limits.** O(k₀ · m) index probes for `m` associations with `k₀` keys in
the first. A malformed argument (not a list of associations or rules) yields
`KeyIntersection::invar` and leaves the call unevaluated. `KeyComplement` is the
sibling that keeps the first association's keys found in *none* of the others.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), ch. 11 (open-addressing hash tables).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

Each association is cut down to the keys shared by all of them, in the key order of
the *first*; the values come from each association unchanged, so the common keys can
differ in value across the result. `KeyComplement` is the complementary operation —
the first association's entries whose keys appear in none of the others.
