### Worked examples

```mathematica
In[1]:= SelectFirst[{1, 3, 5, 6, 7}, EvenQ]  (* the first element passing the predicate *)
```

```mathematica
In[1]:= SelectFirst[Range[10], # > 5 &]  (* a pure function as the predicate *)
```

```mathematica
In[1]:= SelectFirst[{1, 3, 5}, EvenQ, None]  (* no match, so the supplied default *)
```

```mathematica
In[1]:= SelectFirst[{1, 3, 5}, EvenQ]  (* no match and no default: Missing["NotFound"] *)
```

### Notes

`SelectFirst[list, pred]` returns the first element `e` with `pred[e]` equal to
`True`, scanning left to right and stopping at the first hit — the single-element
companion to `Select`. With no match it returns `Missing["NotFound"]`, or the
`default` given as a third argument. Over an association the values are tested and
the matching value is returned.
