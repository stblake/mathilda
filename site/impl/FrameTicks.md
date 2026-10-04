---
source: src/graphics/render.c
---
**Definition.** `FrameTicks` is an option for `Graphics` and `Plot` that specifies
the tick marks on the edges of a frame. `FrameTicks -> Automatic` (the default)
draws major and minor ticks with labels on every drawn frame edge;
`FrameTicks -> None` keeps the frame box but draws no ticks; the
`{{left, right}, {bottom, top}}` form selects `Automatic` or `None` per edge. It is
an inert option keyword — no builtin and no value of its own (it carries no
attributes), declared with its docstring in `info.c`.

**Representation.** A bare `EXPR_SYMBOL`. It is read as the left of a
`FrameTicks -> value` rule by the 2D renderer in `src/graphics/render.c`: the option
loop matches `name == SYM_FrameTicks` and interprets the value — `Automatic`/`True`
turns ticks on for every drawn edge, `None`/`False` keeps the frame box but
suppresses ticks and labels, and a `{{l, r}, {b, t}}` nesting is read per edge.

**Usage & limits.** Meaningful only inside a `Graphics`/`Plot` call that draws a
frame; with `Frame -> False` there are no frame edges to tick. The symbol itself
performs no computation — all behaviour is in the renderer's reading of the rule.
