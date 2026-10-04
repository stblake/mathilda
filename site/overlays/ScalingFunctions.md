### Worked examples

```mathematica
In[1]:= ScalingFunctions -> {"Log", None}  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= MemberQ[Attributes[ScalingFunctions], Protected]  (* protected, with no value of its own *)
```

```mathematica
In[1]:= Head[Plot[Exp[x], {x, 1, 10}, ScalingFunctions -> "Log"]]  (* accepted by its owning plot builtin *)
```

### Notes

`ScalingFunctions` is an inert option keyword for `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot` and `StreamPlot` — a coordinate transform applied to one or both axes ("Log", "Log10", "Log2", "Reverse", or a per-axis pair). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ScalingFunctions -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
