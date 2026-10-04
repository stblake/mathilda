---
source: src/graphics/contourplot.c
---
**Definition.** `ContourLabels` is an **option symbol** for `ContourPlot` that controls
whether each contour level is annotated with its `z`-value. It has no builtin and no
value — an inert, `Protected` keyword read out of `ContourPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_ContourLabels`). `ContourPlot`'s
option handler (`src/graphics/contourplot.c`) reads a `ContourLabels -> b` rule: `True`
draws `z`-value text at the first visible point of each contour level, and the default
`False` draws no labels.

**Usage & limits.** `Protected`. Meaningful only inside a `ContourPlot` call; on its own
it evaluates to itself, and it is read directly from the option sequence rather than
from an `Options[ContourPlot]` default. Related keywords are `Contours`, `ContourStyle`
and `ContourShading`.
