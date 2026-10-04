### Worked examples

```mathematica
In[1]:= Insert[{a, b, c}, x, 2]  (* x goes before position 2 *)
```

```mathematica
In[1]:= Insert[{a, b, c}, x, -1]  (* -1 appends at the end *)
```

```mathematica
In[1]:= Insert[{a, b, c, d}, x, {{2}, {4}}]  (* insert at several positions at once *)
```

### Notes

`Insert[expr, elem, pos]` inserts `elem` so that it occupies position `pos` in
the result, shifting later elements along; negative positions count from the end,
so `pos = -1` appends. A position may be a path for nested insertion, and a list
of positions `{{p1}, {p2}, ...}` inserts a copy of `elem` at each (positions
refer to the original expression). The original head is preserved, so `Insert`
works on any expression, not just lists.
