### Worked examples

```mathematica
In[1]:= KeySort[<|c -> 3, a -> 1, b -> 2|>]  (* canonical order of the keys *)
In[2]:= KeySort[<|3 -> x, 1 -> y, 2 -> z|>, Greater]  (* a custom ordering function *)
```

### Notes

The one-argument form uses Mathilda's canonical order (`expr_compare`); because keys
are distinct it is a total order. `KeySort[assoc, p]` reuses `Sort`/`Ordering` with
the ordering function `p`, so an ordering `p` cannot decide — such as `Greater` on
symbolic keys — leaves those keys in their original order. Use `KeySortBy` to sort by
a function of each key.
