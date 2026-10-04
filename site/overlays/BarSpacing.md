### Worked examples

```mathematica
In[1]:= BarSpacing -> 0.3  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[BarSpacing], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[BarChart[{3, 1, 2}, BarSpacing -> 0.5]]  (* accepted by its owning plot builtin *)
```

### Notes

`BarSpacing` is an inert option keyword for `BarChart` and `Histogram` — the gap between bars as a fraction of bar width (default 0.2; 0 gives touching bars, 1 all gap). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `BarSpacing -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
