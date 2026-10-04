# Association

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Association[key1 -> val1, key2 -> val2, ...]  (also written <|...|>)`**

Represents an association mapping keys to values with unique, insertion-ordered keys (last value wins on duplicates). Arguments may be rules, lists of rules, or other associations.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= <|"a" -> 1, "b" -> 2|>
Out[1]= <|"a" -> 1, "b" -> 2|>

In[2]:= <|"a" -> 1, "b" -> 2, "a" -> 99|>
Out[2]= <|"a" -> 99, "b" -> 2|>

In[3]:= <|"a" -> 10, "b" -> 20|>[["b"]]
Out[3]= 20

In[4]:= <|"a" -> 10, "b" -> 20|>["a"]
Out[4]= 10

In[5]:= <|"a" -> <|"b" -> 5|>|>["a", "b"]
Out[5]= 5

In[6]:= Association[{{"a" -> 1}, {"b" -> 2, {"c" -> 3}}}]
Out[6]= <|"a" -> 1, "b" -> 2, "c" -> 3|>
```

### Applications (3)

```mathematica
In[7]:= Association[a -> 1, b -> 2]
Out[7]= <|a -> 1, b -> 2|>
```

Duplicate key: last value wins, first position kept

```mathematica
In[8]:= <|a -> 1, b -> 2, a -> 3|>
Out[8]= <|a -> 3, b -> 2|>
```

Lists of rules and associations are spliced

```mathematica
In[9]:= Association[{a -> 1, b -> 2}, <|c -> 3|>]
Out[9]= <|a -> 1, b -> 2, c -> 3|>
```

## Implementation notes

**Algorithm.** `builtin_association` normalises its arguments — any mix of
`Rule`/`RuleDelayed` nodes, `List`s of such rules, and existing associations
(which are spliced) — into one flat rule array via `collect_entries`, then hands
it to `assoc_from_rules`. That pass de-duplicates keys with the rule that *first
occurrence fixes position, last occurrence fixes value*, keeping each entry's own
head so a `RuleDelayed` stays delayed. When the input is already canonical (all
direct 2-arg rules, no splicing, no duplicate keys) the builtin returns `NULL`,
so the common `<|a -> 1, b -> 2|>` literal costs nothing to re-evaluate.

**Data structures.** De-duplication is driven by a transient `KeyIndex`: a
fixed-capacity open-addressing hash set over borrowed `Expr*` keys, sized once to
the maximum key count (power of two, load factor < 0.5) so it never rehashes and
every probe is branch-predictable. Keys are compared with `expr_eq` and hashed
with `expr_hash`. A separate persistent `AssocIndex` (`assoc_index.c`) is
attached lazily to the surviving node on the first single-key read, not at
construction — the fixed-point evaluator keeps the original literal node and
discards the rebuilt one, so an eagerly-attached index would be thrown away.

**Complexity / limits.** Construction and de-duplication are amortised `O(n)` in
the number of entries. Once the lazy index exists, `Lookup`/`KeyExistsQ` and
`assoc[key]` are `O(1)` amortised; without it a single-key probe falls back to an
`O(n)` linear `assoc_scan`. A malformed argument (anything not a rule / rule-list
/ association) leaves the node unevaluated.

**Attributes:** `Protected`.

## References

- T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing).
- Source: [`src/assoc.c`](https://github.com/stblake/mathilda/blob/main/src/assoc.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)

## Notes & additional examples

### Notes

An `Association` maps keys to values with unique, insertion-ordered keys and is
written `<|k1 -> v1, ...|>`. Its arguments may be rules, lists of rules, or other
associations, all flattened into one canonical set; on a duplicate key the first
occurrence fixes the position and the last fixes the value. Lookups by key are
`O(1)` amortised — the association carries a hash index built on first use — so
associations are the right structure for keyed records and counters. `Keys`,
`Values`, `Lookup`, `KeyDrop`, and `Merge` all operate on them.
