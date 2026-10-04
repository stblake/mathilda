### Worked examples

```mathematica
In[1]:= Key["a"][<|"a" -> 1, "b" -> 2|>]  (* the curried lookup operator *)
```

```mathematica
In[1]:= <|a -> 10, b -> 20|>[[Key[a]]]  (* Key as a Part specification *)
```

```mathematica
In[1]:= SortBy[{<|"n" -> 3|>, <|"n" -> 1|>, <|"n" -> 2|>}, Key["n"]]  (* name a field in a record pipeline *)
```

### Notes

`Key[k]` names the key `k` of an association and is inert on its own. It acquires
meaning in three places: the operator form `Key[k][assoc]` returns the value at
`k` (or `Missing["KeyAbsent", k]`); as a part specification `assoc[[Key[k]]]` it
reaches that entry; and anywhere a key is read — `Lookup`, `KeyDrop`,
`JoinAcross` — `Key[k]` unwraps to `k`. Its purpose is to let a key act as an
ordinary operator, so record pipelines such as `SortBy[rows, Key["field"]]` and
`GroupBy[rows, Key["field"]]` can select a field by name. Wrapping is needed only
when the bare key would otherwise evaluate to something else.
