---
source: src/graphics/graphics_init.c
---
**Definition.** `LightPink` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightPink", rgb_color(1, 0.925, 0.925))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[1, 0.925, 0.925]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightPink` evaluates to `RGBColor[1, 0.925, 0.925]`, so
`FullForm[LightPink]` is `RGBColor[1, 0.925, 0.925]` and `Head[LightPink]` is `RGBColor`.
The exact `1` red component prints as an integer and the `0.925` green/blue as reals. It
is a distinct constant from `Pink` (`RGBColor[1, 0.5, 0.5]`), one of the pale `Light*`
pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightPink, Disk[]}]`, `PlotStyle -> LightPink`, or a pale fill/background. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.
