### Worked examples

```mathematica
In[1]:= Discard[{1, 2, 3, 4, 5}, EvenQ]  (* the complement of Select *)
```

```mathematica
In[1]:= Discard[Range[10], # > 5 &]
```

```mathematica
In[1]:= Discard[{1, 2, 3, 4, 5}, OddQ, 1]  (* drop at most the first match *)
```

```mathematica
In[1]:= Discard[<|a -> 1, b -> 2, c -> 3|>, OddQ]  (* over an association, tests the values *)
```

### Notes

`Discard[expr, crit]` keeps the elements for which `crit` does *not* return
`True` — exactly the complement of `Select[expr, crit]`. The optional third
argument `Discard[expr, crit, n]` discards at most the first `n` matching
elements and keeps everything after. It works on any non-atomic expression and,
over an association, tests the values while preserving keys.
