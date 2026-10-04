---
source: src/graphics/graphics_init.c
---
**Definition.** `Green` is one of Mathilda's named colour constants, bound as an
**OwnValue** rather than registered as an inert head: `register_color("Green",
rgb_color(0, 1, 0))` in `graphics_init` (`src/graphics/graphics_init.c`) makes the bare
symbol `Green` evaluate to the literal `RGBColor[0, 1, 0]` and marks it `Protected`. The
docstring lives centrally in `src/info.c` (set in `info_init`, which runs first).

**Representation.** `Green` evaluates to `RGBColor[0, 1, 0]`, so `FullForm[Green]` is
`RGBColor[0, 1, 0]` and `Head[Green]` is `RGBColor`. The `color_component` helper emits an
exact `0`/`1` as an integer and anything else as a real, matching Mathematica's
`InputForm`. `Green` is one of the saturated primaries built with `RGBColor` (the
greyscale names use `GrayLevel`).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Green, Line[...]}]`, `PlotStyle -> Green`, `Directive[Green, Dashed]`. `Plot`
evaluates its option values once before storing them, so a held option like
`Epilog -> Green` also sees `RGBColor[0, 1, 0]`. It carries no opacity; use `Opacity[a]` or
the four-argument `RGBColor` for translucency.
