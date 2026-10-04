### Worked examples

```mathematica
In[1]:= StringRiffle[{"a", "b", "c"}]  (* default scheme: a single space *)
```

```mathematica
In[1]:= StringRiffle[{"2024", "01", "02"}, "-"]  (* a custom separator *)
```

```mathematica
In[1]:= StringRiffle[{"x", "y", "z"}, {"(", ", ", ")"}]  (* a {left, sep, right} triple *)
```

### Notes

`StringRiffle` is the inverse of `StringSplit`: it joins a (possibly nested) list
with separators. A plain string separator goes between the top-level elements; a
3-string list is a `{left, sep, right}` delimiter triple.

The default scheme uses a single space at the innermost level and one extra
newline per level above it. Non-string leaves are rendered with `ToString` (so
`27` becomes `"27"`); string leaves are used verbatim.
