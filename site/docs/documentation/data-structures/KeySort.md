# KeySort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KeySort[assoc]`**

Sorts an association into canonical key order.

**`KeySort[assoc, p]`**

Sorts the entries by key using the ordering function p.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= KeySort[<|"c" -> 3, "a" -> 1, "b" -> 2|>]
Out[1]= <|"a" -> 1, "b" -> 2, "c" -> 3|>

In[2]:= KeySort[<|"c" -> 1, "a" -> 2, "b" -> 3|>, -Order[#1, #2] &]
Out[2]= <|"c" -> 1, "b" -> 3, "a" -> 2|>
```

### Applications (2)

Canonical order of the keys

```mathematica
In[3]:= KeySort[<|c -> 3, a -> 1, b -> 2|>]
Out[3]= <|a -> 1, b -> 2, c -> 3|>
```

A custom ordering function

```mathematica
In[4]:= KeySort[<|3 -> x, 1 -> y, 2 -> z|>, Greater]
Out[4]= <|3 -> x, 2 -> z, 1 -> y|>
```

## Implementation notes

**Algorithm.** `builtin_keysort` with one argument copies the entries and `qsort`s
them with `rule_key_cmp`, which orders by `expr_compare` of the keys. Keys are
distinct, so this is a total order and the sort is well defined. The two-argument
`KeySort[assoc, p]` orders by a user ordering function `p`: it extracts the keys and
calls `Ordering[keys, All, p]`, then reindexes the entries by that permutation —
sharing `Sort`'s merge sort, so ties and a non-boolean `p` behave exactly as `Sort`
would.

**Data structures.** An `Expr**` array of entry copies; the system `qsort` (1-arg) or
the `Ordering` permutation vector (2-arg) rebuilt into a fresh `Association`.

**Complexity / limits.** O(n log n). When `p` is undecidable on the keys (e.g.
`Greater` on bare symbols) `Ordering` leaves them in place, so the result keeps the
original order rather than erroring.

- `KeySort[assoc, p]` orders the keys exactly as `Sort[Keys[assoc], p]` would
  (the shared merge sort), so a `p` that does not evaluate to `False`/`-1`
  leaves the order alone, as in Mathematica 15.

**Attributes:** `Protected`.

## References

- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

The one-argument form uses Mathilda's canonical order (`expr_compare`); because keys
are distinct it is a total order. `KeySort[assoc, p]` reuses `Sort`/`Ordering` with
the ordering function `p`, so an ordering `p` cannot decide — such as `Greater` on
symbolic keys — leaves those keys in their original order. Use `KeySortBy` to sort by
a function of each key.
