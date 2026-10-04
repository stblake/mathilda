---
source: src/graphics/plot.c
---
**Definition.** `$PlotResample` is an **internal metadata head** that `Plot` embeds inside
the `Graphics` object it returns, so the renderer can re-sample the curves when the view is
zoomed. It is not a user function: `graphics_init` (`src/graphics/graphics_init.c`) stamps
it `HoldAll | Protected` and sets a docstring that says outright it is not intended for
direct use; the node is built by `make_resample_meta` in `src/graphics/plot.c` and consumed
by the renderer in `src/graphics/show.c`.

**Representation.** `$PlotResample` is a bare `EXPR_SYMBOL` (`Head[$PlotResample]` is
`Symbol`, `Attributes[$PlotResample]` is `{HoldAll, Protected}`). In use it heads a
three-argument form `$PlotResample[var, {bodies}, {opts...}]`, where `bodies` are the
unevaluated plotted expressions and `opts` is the packed option vector (`plotPoints`,
`maxRecursion`, `maxPlotPoints`, mesh, region function, exclusions, colour function and its
scaling, filling, filling style). `HoldAll` is what keeps the bodies and the held option
values unevaluated through the surrounding `Graphics` re-evaluation.

**Usage & limits.** It exists only as a carrier: `show.c` recognises the head while drawing
and, on a zoom, re-samples the stored bodies over the new coordinate window instead of
rescaling a fixed polyline. It has no standalone value — written by hand it simply stays
unevaluated — and is an implementation detail of `Plot`'s zoomable output.
