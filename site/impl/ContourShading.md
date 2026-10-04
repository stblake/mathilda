---
source: src/graphics/graphics_init.c
---
**Definition.** `ContourShading` is an **option name** for `ContourPlot` that controls
whether the grid cells between contour lines are filled. It is a bare option symbol, not a
function: `graphics_init` (`src/graphics/graphics_init.c`) stamps it `Protected` and sets
its docstring, but it has no builtin and no DownValues. The `ContourPlot` sampler reads it
from the call's option list.

**Representation.** `ContourShading` stays an inert, `Protected` `EXPR_SYMBOL`, appearing
only on the left of a rule; `FullForm[ContourShading -> False]` is
`Rule[ContourShading, False]`. The value is `True`, `False`, or `Automatic`.

**Usage & limits.** With `ContourShading -> True` each cell is filled by its `z` value,
using the supplied `ColorFunction` or the built-in blue-cyan-yellow-red thermal ramp;
`False` draws only the contour lines; `Automatic` enables shading only when a
`ColorFunction` is set. It is consumed at render time and affects only the fill, never the
sampled data. It is specific to `ContourPlot`, not a general `Graphics` directive.
