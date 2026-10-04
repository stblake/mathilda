### Worked examples

```mathematica
In[1]:= JoinAcross[{<|"id" -> 1, "x" -> a|>, <|"id" -> 2, "x" -> b|>}, {<|"id" -> 1, "y" -> p|>, <|"id" -> 1, "y" -> q|>}, "id"]
```

```mathematica
In[1]:= JoinAcross[{<|"id" -> 1, "x" -> a|>, <|"id" -> 2, "x" -> b|>}, {<|"id" -> 1, "y" -> p|>}, "id", "Left"]
```

### Notes

`JoinAcross[left, right, spec]` is the relational join over two lists of
associations (rows), matching rows whose join keys agree and merging the matched
pairs into combined records. The join spec is a key, `Key[k]`, a renaming
`k1 -> k2` (different column names on the two sides), or a list of these to join
on several columns. The default is an inner join, so unmatched rows are dropped;
a fourth argument `"Left"`, `"Right"`, or `"Outer"` keeps them, filling the
absent columns with `Missing["Unmatched"]`. The right rows are hash-indexed on
the join key, so the join is linear in the combined size.
