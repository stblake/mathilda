### Worked examples

```mathematica
In[1]:= MinimalBy[{5, -2, 3, 2}, Abs]  (* both -2 and 2 tie for the smallest absolute value *)
```

```mathematica
In[1]:= MinimalBy[{1, -5, 3, -2}, Abs]  (* a unique minimum is still returned as a list *)
```

```mathematica
In[1]:= MinimalBy[Abs][{5, -2, 3, 2}]  (* operator form MinimalBy[f][list] *)
```

### Notes

`MinimalBy[list, f]` returns the element (or elements) of `list` for which `f` is
smallest, by Mathilda's canonical order — the dual of `MaximalBy`. All ties are
returned, in original order, so the result is always a list (`{-2, 2}` both have
absolute value 2). The operator form `MinimalBy[f]` defers the collection; over an
association the entries whose *value* minimises `f` are returned.
