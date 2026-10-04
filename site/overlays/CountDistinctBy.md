### Worked examples

```mathematica
In[1]:= CountDistinctBy[{1, 2, 3, 4, 5, 6}, EvenQ]  (* two classes: even and odd *)
```

```mathematica
In[1]:= CountDistinctBy[{-1, 1, -2, 3}, Abs]  (* |x| takes three distinct values *)
```

```mathematica
In[1]:= CountDistinctBy[Range[10], Mod[#, 3] &]
```

### Notes

`CountDistinctBy[list, f]` counts the distinct values of `f[element]` — the number
of groups `GatherBy[list, f]` would produce, without materialising the groups. It
is a one-pass hash count: `f` is applied once per element and the distinct results
are tallied. Use it to ask how many categories a key function induces.
