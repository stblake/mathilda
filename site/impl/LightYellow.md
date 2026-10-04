---
source: src/graphics/graphics_init.c
---
**Definition.** `LightYellow` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightYellow", rgb_color(1, 1, 0.85))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[1, 1, 0.85]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightYellow` evaluates to `RGBColor[1, 1, 0.85]`, so
`FullForm[LightYellow]` is `RGBColor[1, 1, 0.85]` and `Head[LightYellow]` is `RGBColor`.
The exact `1` red/green components print as integers and the `0.85` blue as a real
(`color_component`). It is a distinct constant from `Yellow` (`RGBColor[1, 1, 0]`), one of
the pale `Light*` pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightYellow, Disk[]}]`, `PlotStyle -> LightYellow`, a pale highlight fill. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.
