### Worked examples

```mathematica
In[1]:= VectorPoints -> 20  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[VectorPoints], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 5]]  (* accepted by its owning plot builtin *)
```

### Notes

`VectorPoints` is an inert option keyword for `VectorPlot` — an integer n giving an n*n seed grid of arrows (default 15). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `VectorPoints -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
