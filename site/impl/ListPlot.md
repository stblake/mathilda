---
source: src/graphics/listplot.c
---
**Algorithm.** `builtin_listplot` is a plain `Protected` builtin — **not**
`HoldAll` — so its data is already evaluated (`ListPlot[Range[10]]`,
`ListPlot[Table[i^2, {i, 5}]]` and named colours in options all work). `classify`
/ `extract_dataset` sort the argument into one or more datasets: a flat list with
no sublists is read as **heights** `{y_i}` plotted at `x` running over `DataRange`
(default `1..n`), with non-numeric entries skipped but still consuming their index
slot; a list in which *every* element is a numeric 2-pair is read as **explicit
points** `{{x, y}, ...}`; anything else (some elements are sublists) is read as
**multiple datasets**, one per sublist, in distinct `palette_color`s.
`DataRange -> All` forces the multi reading of a pair list. Numeric coercion goes
through `to_double` (literal `Integer`/`Real` fast path, else `N[]`). Each dataset
is emitted by `emit_dataset` as a coordinate `List` wrapped in `Line[...]` when
`Joined -> True` (default `False`) and otherwise `Point[...]`; multi-dataset
styling gets a per-dataset palette colour with non-colour directives in their own
`List` scope. `Filling` draws continuous quads under a joined curve
(`gfx_build_fill_quads`) or vertical stems under unjoined points;
`ScalingFunctions` transform world coordinates in place; `PlotLegends` attaches
legend metadata. The result is an inert `Graphics[prims, opts...]`.

**Data structures.** `PointSet {double* xs; double* ys; size_t n}`, one per live
dataset (`classify` returns an array of borrowed dataset `Expr*`); the option
bundle from `split_options`.

**Complexity / limits.** `O(total points)` — there is no sampling, the data is
plotted as given. Default `Axes -> True`, `AspectRatio -> 1/GoldenRatio`
(distinct from the field plotters' `1`), `PlotStyle` colour
`RGBColor[0.2, 0.4, 0.8]`. The default marker size is resolved by the renderer
from scatter density rather than emitted here; `PlotMarkers` is accepted but its
glyphs are not yet drawn.
