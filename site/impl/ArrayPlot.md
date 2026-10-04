---
source: src/graphics/arrayplot.c
---
**Algorithm.** `builtin_arrayplot` is a plain `Protected` builtin — **not**
`HoldAll` — so its argument is an already-evaluated array, rendered as a discrete
grid of coloured cells with **no interpolation** (unlike `DensityPlot`). Two
ingestion paths: an `EXPR_NDArray` is loaded by `na_load_matrix`
(`src/linalg/numarray.h`, a dense numeric buffer, every cell numeric); a nested
`List` goes to `arrayplot_load`, which walks the rectangular list-of-lists and
classifies each cell as either a **colour literal**
(`RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`, painted directly) or a **plain number**
(`expr_to_real_double`). Numeric and colour cells may mix freely; a ragged shape,
a non-list-of-lists, or a cell that is neither colour nor number makes the whole
call decline (`NULL`). Each numeric cell is coloured by `ap_color`: the default
`ColorFunction` is **"Greyscale"** (white at the array minimum, black at the
maximum — matching Mathematica), with `ColorFunctionScaling -> True` normalising
`(v - zmin)/zspan`; `ColorRules -> {key -> colour}` overrides exact-matching
values before the ramp. Row 0 is drawn at the **top** (`y0 = rows-1-i`), column 0
at the **left**; each cell is a colour directive plus a `Rectangle`, overlapping
its later-drawn neighbours to close seams. `Mesh -> All` adds `rows+1` horizontal
and `cols+1` vertical `GrayLevel[0.5]` grid `Line`s. `PlotLegends -> Automatic`
appends a `$StreamColorBar` only when at least one cell is numeric. The result is
an inert `Graphics[prims, opts...]`.

**Data structures.** `ArrayPlotOpts`; parallel arrays `const Expr** cell_colors`
(borrowed colour-literal pointers, or NULL) and `double* buf` (values), each
`rows*cols`; `Expr** prims` sized `rows*cols*2 + mesh + 4`. `AspectRatio`
defaults to `rows/cols` (square cells); `PlotRange` is the exact extent
`{{0, cols}, {0, rows}}`.

**Complexity / limits.** `O(rows*cols)`, no sampling — a literal heatmap or, when
every cell is a colour, a raw pixel-grid renderer. Defaults `Frame -> True`,
`Axes -> False`.
