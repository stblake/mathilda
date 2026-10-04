# KeyDrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeyDrop[assoc, key]  |  KeyDrop[assoc, {k1, ...}]`**

Gives assoc with the specified keys removed (order preserved).

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, "b"]
Out[1]= <|"a" -> 1, "c" -> 3|>

In[2]:= KeyDrop[<|"a" -> 1, "b" -> 2, "c" -> 3|>, {"a", "c"}]
Out[2]= <|"b" -> 2|>
```

### Scope (1)

```mathematica
In[3]:= KeyDrop[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 3, "b" -> 4|>}, "b"]
Out[3]= {<|"a" -> 1|>, <|"a" -> 3|>}
```

### Applications (3)

```mathematica
In[4]:= KeyDrop[<|a -> 1, b -> 2, c -> 3|>, b]
Out[4]= <|a -> 1, c -> 3|>
```

Drop several keys at once

```mathematica
In[5]:= KeyDrop[<|a -> 1, b -> 2, c -> 3|>, {a, c}]
Out[5]= <|b -> 2|>
```

Threads over a list of records

```mathematica
In[6]:= KeyDrop[{<|a -> 1, b -> 2|>, <|a -> 3, b -> 4|>}, a]
Out[6]= {<|b -> 2|>, <|b -> 4|>}
```

## Implementation notes

**Algorithm.** `builtin_keydrop` is the drop half of `key_drop_take`, the shared
core with `KeyTake`. Given a key or a `List` of keys it builds a `KeyIndex` of
that drop set once, then rebuilds the association keeping exactly the entries
whose key is *absent* from the set — order preserved. A `List` of associations is
threaded element-wise (the column-of-records form), each element delegated back
to the same builtin.

**Data structures.** A transient `KeyIndex` open-addressing hash set over the
requested keys; the surviving entries are deep-copied into a fresh canonical
`Association`. The compiled evaluator calls the lower-level `assoc_key_select`
directly, with no call-node round-trip.

**Complexity / limits.** `O(n + m)` for `n` entries and `m` requested keys,
versus the `O(n·m)` of repeated membership scans. Returns `NULL` unless the first
argument is an association (or a non-empty list of them).

**Attributes:** `Protected`.

## References

**See also:** [KeyTake](../../data-structures/KeyTake/)

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile_assoc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_assoc.c)

## Notes & additional examples

### Notes

`KeyDrop[assoc, key]` returns the association with the given key removed, keeping
the order of the remaining entries; the second argument may be a single key or a
list of keys. Given a list of associations it threads — dropping the key(s) from
each record, which is the column-of-records form. It is the complement of
`KeyTake`, and both build a hash index of the key set once so the work is linear
rather than quadratic in the key count.
