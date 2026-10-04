### Worked examples

```mathematica
In[1]:= Catenate[{{1, 2}, {3, 4}, {5}}]
```

```mathematica
In[1]:= Catenate[{{1}, {2, 3}, {4, 5, 6}}]  (* ragged lengths are fine *)
```

```mathematica
In[1]:= Catenate[{<|a -> 1|>, <|b -> 2|>}]  (* associations take part through their values *)
```

```mathematica
In[1]:= Catenate[<|x -> {1, 2}, y -> {3}|>]  (* the collection may itself be an association *)
```

### Notes

`Catenate[{e1, e2, ...}]` flattens one level: the parts must share a head and
their elements are concatenated under it, so `Catenate[{{1,2},{3,4}}]` is
`{1,2,3,4}`. It differs from `Join`, which concatenates several *arguments*
rather than one list of parts. Associations join through their values
(`Catenate[{<|a->1|>, <|b->2|>}]` is `{1, 2}`), and the outer collection may
itself be an association. A packed matrix catenates as a reshape of its row-major
buffer, so the operation stays on the buffer.
