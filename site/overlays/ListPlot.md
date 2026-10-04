### Worked examples

```mathematica
(* a flat list of heights, plotted at {i, y_i} *)
In[1]:= ListPlot[{1, 4, 9, 16, 25}]
```

```mathematica
(* a list of coordinate pairs is a scatter plot *)
In[1]:= ListPlot[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}]
```

```mathematica
(* Joined -> True connects the points with one Line *)
In[1]:= ListPlot[Table[Sin[n], {n, 20}], Joined -> True]
```

```mathematica
(* unlike Plot, ListPlot is not HoldAll, so Range/Table are evaluated first *)
In[1]:= Length[Cases[ListPlot[Range[10]], _Point, Infinity]]
```

```mathematica
(* the default aspect ratio, distinct from the field plotters' 1 *)
In[1]:= Options[ListPlot, AspectRatio]
```

### Notes

`ListPlot` is `Protected` but **not** `HoldAll`, so its data is evaluated before
classification — `ListPlot[Range[10]]` and `ListPlot[Table[i^2, {i, 5}]]` both
work. The argument is read as heights (a flat list, plotted at `{i, y_i}` over
`DataRange`), explicit coordinate pairs (every element a numeric 2-list), or
several datasets (some elements are sublists), each in a distinct `ColorData[97]`
palette colour. Non-numeric entries are skipped but still consume their index
slot.

`Joined -> True` emits one `Line` per dataset instead of a `Point` cloud;
`Filling` draws continuous quads under a joined curve or vertical stems under
points. The default `AspectRatio` is `1/GoldenRatio` (the field plotters default
to `1`), and the point marker size is chosen by the renderer from scatter density
unless an explicit `PointSize` overrides it.
