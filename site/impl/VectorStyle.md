---
source: src/graphics/graphics_init.c
---
**Definition.** `VectorStyle` is an **option name** for `VectorPlot` giving the style
directive(s) applied globally to every arrow. It is a bare option symbol, not a function:
`graphics_init` (`src/graphics/graphics_init.c`) stamps it `Protected` and sets its
docstring, but it has no builtin and no DownValues. The `VectorPlot` renderer reads it from
the call's option list.

**Representation.** `VectorStyle` stays an inert, `Protected` `EXPR_SYMBOL`, appearing only
on the left of a rule; its value is one or more style directives (`RGBColor`, `Thickness`,
…). Because named colours are OwnValues, `FullForm[VectorStyle -> Red]` resolves to
`Rule[VectorStyle, RGBColor[1, 0, 0]]`.

**Usage & limits.** Supplied as
`VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorStyle -> Red]`, it applies the given
directive to all arrows and overrides any per-arrow `ColorFunction`. It is consumed at
render time and affects only the arrows' appearance, never the sampled field. It is
specific to `VectorPlot`, not a general `Graphics` directive.
