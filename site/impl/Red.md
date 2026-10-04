---
source: src/graphics/graphics_init.c
---
**Definition.** `Red` is one of Mathilda's named colour constants. It is not an inert
head but an **OwnValue**: `register_color("Red", rgb_color(1, 0, 0))` in `graphics_init`
(`src/graphics/graphics_init.c`) binds the bare symbol `Red` to the literal
`RGBColor[1, 0, 0]` and stamps it `Protected`. The docstring is kept centrally in
`src/info.c`, beside `CMYKColor`, because `info_init` runs before `graphics_init`.

**Representation.** Evaluating `Red` yields `RGBColor[1, 0, 0]`, so `FullForm[Red]` is
`RGBColor[1, 0, 0]` and `Head[Red]` is `RGBColor` — the OwnValue has already fired by the
time anything inspects the result. The helper `color_component` prints an exact `0` or `1`
as a plain integer (matching Mathematica's own `InputForm[Red]`) and any other value as a
real. The saturated primaries use `RGBColor`; the greyscale names (`Black`, `White`,
`Gray`, `LightGray`) use `GrayLevel`.

**Usage & limits.** Because it resolves to a genuine colour literal wherever a bare
argument is evaluated, `Red` works anywhere a style directive is expected —
`Graphics[{Red, Disk[]}]`, `PlotStyle -> Red`, an `ArrayPlot` cell, `Directive[Red, Thick]`
— exactly as the long `RGBColor[...]` form would. `Plot` evaluates each option value once
before storing it (`split_options` in `plot.c`), so `Epilog -> {Red, ...}` also sees the
resolved literal. `Red` carries no opacity; use `Opacity[a]` or `RGBColor[r, g, b, a]` for
a translucent colour.
