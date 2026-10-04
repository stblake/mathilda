### Worked examples

```mathematica
In[1]:= Reap[Scan[Sow, {1, 2, 3}]]  (* Scan visits each element; Sow records the order, Scan itself returns Null *)
```

```mathematica
In[1]:= Reap[Scan[Sow, {{1, 2}, {3, 4}}, {2}]]  (* a level spec reaches the inner elements *)
```

```mathematica
In[1]:= Reap[Scan[Sow, a + b c, {-1}]]  (* depth-first, leaves before roots: the atoms a, b, c *)
```

```mathematica
In[1]:= Scan[Print, {a, b}]  (* applied purely for the side effect; the value is Null *)
```

### Notes

`Scan[f, expr]` applies `f` to each element the way `Map` does but keeps nothing:
it is for side effects (printing, sowing, logging) and always returns `Null`. The
traversal is depth-first with **leaves before roots**, so a level spec such as
`{-1}` visits the atoms of an expression before the subexpressions that contain
them.

A `Throw` inside `f` exits to an enclosing `Catch`; a `Return[ret]` evaluated
directly as `f` makes `Scan` itself return `ret`. Over an association, `f` is
applied to each value.
