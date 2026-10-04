---
source: src/graphics/graphics_init.c
---
**Definition.** `Brown` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("Brown", rgb_color(0.6, 0.4, 0.2))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol `Brown` evaluate to the literal
`RGBColor[0.6, 0.4, 0.2]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `Brown` evaluates to `RGBColor[0.6, 0.4, 0.2]`, so `FullForm[Brown]`
is `RGBColor[0.6, 0.4, 0.2]` and `Head[Brown]` is `RGBColor`. Its components are not `0`
or `1`, so `color_component` renders each as a machine real — unlike the pure primaries,
whose components print as integers. `Brown` is one of the fuller tertiary tones (with
`Orange`, `Pink`, `Purple`).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Brown, Disk[]}]`, `PlotStyle -> Brown`, `Directive[Brown, Thick]`. `Plot`
evaluates its option values once before storing them, so `Epilog -> Brown` also sees the
resolved literal. The light companion `LightBrown` (`RGBColor[0.94, 0.91, 0.88]`) is a
separate constant, not a programmatic lightening of `Brown`. It carries no opacity; use
`Opacity[a]` or the four-argument `RGBColor` for translucency.
