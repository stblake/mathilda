### Worked examples

```mathematica
In[1]:= ContourStyle -> Red  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[ContourStyle], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, ContourStyle -> Red]]  (* accepted by its owning plot builtin *)
```

### Notes

`ContourStyle` is an inert option keyword for `ContourPlot` — the style directive(s) for contour lines; a list cycles through levels, Automatic colours by height, None suppresses lines. It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ContourStyle -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
