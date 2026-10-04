### Worked examples

```mathematica
In[1]:= Association[a -> 1, b -> 2]
```

```mathematica
In[1]:= <|a -> 1, b -> 2, a -> 3|>  (* duplicate key: last value wins, first position kept *)
```

```mathematica
In[1]:= Association[{a -> 1, b -> 2}, <|c -> 3|>]  (* lists of rules and associations are spliced *)
```

### Notes

An `Association` maps keys to values with unique, insertion-ordered keys and is
written `<|k1 -> v1, ...|>`. Its arguments may be rules, lists of rules, or other
associations, all flattened into one canonical set; on a duplicate key the first
occurrence fixes the position and the last fixes the value. Lookups by key are
`O(1)` amortised — the association carries a hash index built on first use — so
associations are the right structure for keyed records and counters. `Keys`,
`Values`, `Lookup`, `KeyDrop`, and `Merge` all operate on them.
