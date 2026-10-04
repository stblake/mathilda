---
source: src/graphics/contourplot.c
---
**Definition.** `Contours` is an **option symbol** for `ContourPlot` that selects the
contour levels. It has no builtin and no value — it is an inert, `Protected` keyword
read out of `ContourPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Contours`). `ContourPlot`'s
option handler (`src/graphics/contourplot.c`) reads a `Contours -> spec` rule where
`spec` is either an integer `n` (draw `n` automatically chosen, evenly spaced levels)
or an explicit list of `z`-values at which to draw contours. The chosen levels drive
the marching-squares sampler and, when shading is on, the colour ramp.

**Usage & limits.** `Protected`. Meaningful only inside a `ContourPlot` call; on its own
it evaluates to itself. It is read directly from the option sequence rather than from a
registered default in `Options[ContourPlot]`. Related keywords are `ContourStyle`,
`ContourLabels` and `ContourShading`.
