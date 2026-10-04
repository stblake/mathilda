---
source: src/graphics/render.c
---
**Definition.** `$StreamColorBar` is an **internal metadata head** — not a user-facing
function. It records the data needed to draw a colour-scale bar beside a plot:
`$StreamColorBar[min, max]` for a density/contour scale, or
`$StreamColorBar[min, max, colourFnOrName]` for a complex-plot phase wheel. It is
emitted by the plot builders and read back by the renderer.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_StreamColorBar`,
`Protected`). `contourplot.c` appends `$StreamColorBar[zmin, zmax, ...]` to the
`Graphics[...]` option list when `PlotLegends -> Automatic`; `complexplot.c` appends
`$StreamColorBar[-Pi, Pi, cfn_or_"Cyclic"]` for its phase bar; `render.c`'s
`find_stream_colorbar` locates the node in the option list and draws the bar.

**Usage & limits.** Internal plumbing: a user never writes `$StreamColorBar` directly,
and it is meaningful only as an item in a graphics option list. Evaluated on its own it
stays symbolic (no evaluation rule). It is documented here only because it is reachable
from the symbol table (`Names["*"]`).

*Internal symbol — shown below only structurally; it has no standalone mathematical
example.*
