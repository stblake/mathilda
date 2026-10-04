### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomChoice[{a, b, c, d}]  (* one element, chosen uniformly *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomChoice[{1, 2, 3, 4, 5, 6}, 10]  (* ten draws, with replacement *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomChoice[{0, 1}, {3, 3}]  (* a 3x3 array of draws *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomChoice[{0.1, 0.1, 0.8} -> {a, b, c}, 12]  (* weighted: c dominates *)
```

### Notes

`RandomChoice` selects **with replacement**, so the same element can appear more
than once and a count larger than the list is fine. For selection *without*
replacement, see `RandomSample`.

The weighted form `RandomChoice[{w1, ...} -> {e1, ...}]` draws each element with
probability proportional to its weight; weights need not sum to `1` and are found
by inverse-CDF with a binary search. A count `n` returns a flat list and `{n1, ...}`
a nested array, each leaf an independent draw.

Every example is seeded with `SeedRandom`, so the choices shown reproduce exactly.
The bare `RandomChoice[list]` returns a single element (a scalar); the `n`-forms
return lists.
