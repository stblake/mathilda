### Worked examples

```mathematica
In[1]:= Attributes[VectorScale]  (* an inert, Protected option keyword *)
```

```mathematica
In[1]:= VectorScale  (* no value of its own: evaluates to itself *)
```

```mathematica
In[1]:= Head[VectorPlot[{1, x}, {x, -1, 1}, {y, -1, 1}, VectorScale -> 0.5]]  (* read inside a VectorPlot *)
```

### Notes

`VectorScale` is a `VectorPlot` option controlling how arrow length encodes the field.
`Automatic` (the default) draws every arrow at equal length, showing direction only;
`None` makes length proportional to the field magnitude; a positive real `f` sets arrow
length to `f ×` the grid spacing.

It is an inert, `Protected` keyword with no value of its own, read straight from the
`VectorPlot` option sequence (it is not in `Options[VectorPlot]`'s default list, so
`SetOptions` on it has no effect). The companion keywords are `VectorPoints` (seed-grid
density) and `VectorStyle` (global style directives).
