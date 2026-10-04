### Worked examples

```mathematica
In[1]:= ReverseSort[{3, 1, 4, 1, 5, 9, 2, 6}]  (* descending order *)
```

```mathematica
In[1]:= ReverseSort[{"banana", "apple", "cherry"}]  (* canonical order works on strings too *)
```

```mathematica
In[1]:= ReverseSort[<|a -> 2, b -> 2, c -> 1|>]  (* descending by value, ties kept in input order *)
```

### Notes

`ReverseSort[list]` sorts into descending order — `Reverse[Sort[list]]` — using
Mathilda's canonical order, so it ranks numbers, strings, and general expressions
alike. Over an association it sorts the entries by value, descending, keeping
equal values in their input order (so `a -> 2, b -> 2` are not swapped); this is
the Mathematica 15 convention, which is subtly *not* `Reverse[Sort[...]]`.
