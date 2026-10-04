---
source: src/graphics/render.c
---
**Definition.** `Frame` is an **option symbol** for `Graphics` and the plotting heads
(`Plot`, `ListPlot`, `ContourPlot`, `DensityPlot`, `VectorPlot`, `ArrayPlot`, ...) that
decides whether a boxed frame with ticks and labels is drawn around the plot. It has no
builtin and no value of its own — it is a keyword read out of the option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Frame`). The renderer and each
plot's option parser read a `Frame -> spec` rule: `True` boxes all four edges;
`False`/`None` draws no frame; `{{left, right}, {bottom, top}}` toggles each edge
individually. In `Plot` a frame takes the place of the default `Axes`. Several heads
(e.g. `ArrayPlot`, `ComplexPlot`) default `Frame -> True`.

**Usage & limits.** It carries a docstring from `info.c` but, unlike the other graphics
option keywords, is **not** marked `Protected` (`Attributes[Frame]` is `{}`), so it can
in principle be reassigned. It is meaningful only inside a graphics/plot call; on its
own it evaluates to itself. Companion keywords are `FrameTicks`, `FrameLabel`, `Axes`
and `AxesLabel`.
