---
source: src/graphics/graphics_init.c
---
**Definition.** `LightGreen` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightGreen", rgb_color(0.88, 1, 0.88))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[0.88, 1, 0.88]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightGreen` evaluates to `RGBColor[0.88, 1, 0.88]`, so
`FullForm[LightGreen]` is `RGBColor[0.88, 1, 0.88]` and `Head[LightGreen]` is `RGBColor`.
The exact `1` green component prints as an integer while the `0.88` reds print as reals
(`color_component`). It is a separate constant from `Green`, one of the pale `Light*`
pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightGreen, Disk[]}]`, `PlotStyle -> LightGreen`, a light `ArrayPlot` cell. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.
