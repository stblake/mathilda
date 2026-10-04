### Worked examples

```mathematica
(* a rotation field, arrows on a 15x15 grid *)
In[1]:= VectorPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
```

```mathematica
(* VectorScale -> None makes arrow length proportional to magnitude *)
In[1]:= VectorPlot[{x, y}, {x, -2, 2}, {y, -2, 2}, VectorScale -> None]
```

```mathematica
(* a RegionFunction masks the grid to the unit disk *)
In[1]:= VectorPlot[{-y, x}, {x, -1.5, 1.5}, {y, -1.5, 1.5}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 1]]
```

```mathematica
(* one Arrow primitive per grid point: a 6x6 grid is 36 arrows *)
In[1]:= Length[Cases[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 6], _Arrow, Infinity]]
```

```mathematica
(* HoldAll and Protected *)
In[1]:= Attributes[VectorPlot]
```

### Notes

`VectorPlot` is `HoldAll`. It samples `{vx, vy}` on an `N x N` grid (`VectorPoints`
default 15) and draws one `Arrow` per grid point, **centred** on that point. Arrow
sizing is screen-normalised so arrows stay legible across mixed-scale axes:
`VectorScale -> Automatic` (default) gives a fixed length, `None` makes length
proportional to magnitude, and a real value scales relative to the grid spacing.

The default `ColorFunction` is the Viridis ramp keyed to speed; a custom function
is tried as `f[vx, vy, speed]` then `f[speed]`. When a `ScalingFunctions` transform
is active the arrow direction is corrected by a finite-difference Jacobian. The grid
is always square (there is no per-axis point count).
