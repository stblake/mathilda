---
source: src/graphics/render.c
---
**Definition.** `$PlotLegendData` is an **internal metadata head** — not a user-facing
function. It carries the legend of a plot inside the option list of a `Graphics[...]`
object: `$PlotLegendData[{colour, label}, {colour, label}, ...]`, one entry per curve.
It is built by `Plot` (and friends) when `PlotLegends` is given, and read back by the
renderer when it draws the legend box.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_PlotLegendData`) with the
attributes `{HoldAll, Protected}`. `plot.c`'s `build_legend_meta` constructs it from the
`PlotLegends` option and appends it to the `Graphics[...]` option sequence;
`render.c`'s `find_legend_data` locates it there, and `Show[]` concatenates the
`$PlotLegendData` entries of several graphics into one combined legend (while dropping
`$PlotResample`). `HoldAll` keeps the stored `{colour, label}` pairs from re-evaluating.

**Usage & limits.** Internal plumbing: a user never writes `$PlotLegendData` directly,
and it is meaningful only as an item in a graphics option list. Evaluated on its own it
stays symbolic (it has no evaluation rule). It is documented here for completeness
because it is reachable from the symbol table (`Names["*"]`).

*Internal symbol — shown below only structurally; it has no standalone mathematical
example.*
