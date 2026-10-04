# JoinAcross

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`JoinAcross[{a1, ...}, {b1, ...}, spec]`**

Joins two lists of associations, merging each ai with every bj whose join keys agree (inner join). spec is a key, Key\[k\], k1 -\> k2 (keys named differently on the two sides), or a list of these.

**`JoinAcross[{a1, ...}, {b1, ...}, spec, type]`**

type is "Inner", "Left", "Right" or "Outer"; unmatched rows get Missing\["Unmatched"\] for the other side's keys. Option KeyCollisionFunction -\> Left | Right | f resolves a non-join key present on both sides (f\[k\] gives the pair of new keys).

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, Key["id"]]
Out[1]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>}

In[2]:= JoinAcross[{<|"id" -> 1, "name" -> "Ada"|>, <|"id" -> 2, "name" -> "Bob"|>}, {<|"id" -> 1, "dept" -> "R&D"|>, <|"id" -> 3, "dept" -> "Ops"|>}, "id", "Outer"]
Out[2]= {<|"id" -> 1, "name" -> "Ada", "dept" -> "R&D"|>, <|"id" -> 2, "name" -> "Bob", "dept" -> Missing["Unmatched"]|>, <|"id" -> 3, "name" -> Missing["Unmatched"], "dept" -> "Ops"|>}

In[3]:= JoinAcross[{<|"k" -> 1, "x" -> 10|>}, {<|"key" -> 1, "y" -> 20|>}, "k" -> "key"]
Out[3]= {<|"k" -> 1, "x" -> 10, "key" -> 1, "y" -> 20|>}
```

### Options (1)

```mathematica
In[4]:= JoinAcross[{<|"a" -> 1, "v" -> "L"|>}, {<|"a" -> 1, "v" -> "R"|>}, "a", KeyCollisionFunction -> Right]
Out[4]= {<|"a" -> 1, "v" -> "R"|>}
```

### Applications (2)

```mathematica
In[5]:= JoinAcross[{<|"id" -> 1, "x" -> a|>, <|"id" -> 2, "x" -> b|>}, {<|"id" -> 1, "y" -> p|>, <|"id" -> 1, "y" -> q|>}, "id"]
Out[5]= {<|"id" -> 1, "x" -> a, "y" -> p|>, <|"id" -> 1, "x" -> a, "y" -> q|>}

In[6]:= JoinAcross[{<|"id" -> 1, "x" -> a|>, <|"id" -> 2, "x" -> b|>}, {<|"id" -> 1, "y" -> p|>}, "id", "Left"]
Out[6]= {<|"id" -> 1, "x" -> a, "y" -> p|>, <|"id" -> 2, "x" -> b, "y" -> Missing["Unmatched"]|>}
```

## Implementation notes

**Algorithm.** `builtin_joinacross` is the relational join over two lists of
associations (rows). The join spec is a key `k`, `Key[k]`, a renaming
`k1 -> k2`, or a `List` of these (join on several columns at once). The right
rows are hashed on their join-key tuple; each left row's tuple is then probed,
and every matching right row is merged into the left row. A fourth argument
selects the join type (`"Inner"` default, `"Left"`, `"Right"`, `"Outer"`), and
the option `KeyCollisionFunction -> Left | Right | f` decides how colliding
non-join keys are kept.

**Data structures.** A `KeyIndex`/hash table over the right rows keyed by the
owned join-key tuple (`join_tuple` builds a `List` of the key values, returning
`NULL` when a join key is absent so the row never matches). Output order follows
Mathematica 15: matched pairs left-major, then unmatched left rows, then
unmatched right rows; unmatched rows carry the other side's common keys as
`Missing["Unmatched"]`, and disagreeing key sets are padded to the union with
`Missing["NotAvailable"]`.

**Complexity / limits.** `O(|left| + |right| + |output|)` — the hash index
replaces the quadratic scan of nested-loop matching. Invalid arguments route
their diagnostics through `ops_msg` (so `Quiet[]`/`Check[]` are honoured).

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Garcia-Molina, Ullman and Widom, *Database Systems: The Complete Book*, 2nd ed. (Pearson, 2008), §15.4 (hash join).
- Source: [`src/assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/src/assoc_ops.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)

## Notes & additional examples

### Notes

`JoinAcross[left, right, spec]` is the relational join over two lists of
associations (rows), matching rows whose join keys agree and merging the matched
pairs into combined records. The join spec is a key, `Key[k]`, a renaming
`k1 -> k2` (different column names on the two sides), or a list of these to join
on several columns. The default is an inner join, so unmatched rows are dropped;
a fourth argument `"Left"`, `"Right"`, or `"Outer"` keeps them, filling the
absent columns with `Missing["Unmatched"]`. The right rows are hash-indexed on
the join key, so the join is linear in the combined size.
