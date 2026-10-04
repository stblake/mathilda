---
source: src/graphics/graphics_init.c
---
**Definition.** `Cyan` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("Cyan", rgb_color(0, 1, 1))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol `Cyan` evaluate to the literal
`RGBColor[0, 1, 1]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`, beside `CMYKColor`.

**Representation.** `Cyan` evaluates to `RGBColor[0, 1, 1]`, so `FullForm[Cyan]` is
`RGBColor[0, 1, 1]` and `Head[Cyan]` is `RGBColor`. The `color_component` helper prints the
exact `0`/`1` components as integers. `Cyan`, `Magenta` (`RGBColor[1, 0, 1]`) and `Yellow`
(`RGBColor[1, 1, 0]`) are the secondary (subtractive) primaries — the same colour a
`CMYKColor` single-ink directive converts to.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Cyan, Disk[]}]`, `PlotStyle -> Cyan`, `Directive[Cyan, Thick]`. `Plot` evaluates
its option values once before storing them, so `Epilog -> Cyan` also sees
`RGBColor[0, 1, 1]`. It carries no opacity; use `Opacity[a]` or the four-argument
`RGBColor` for translucency.
