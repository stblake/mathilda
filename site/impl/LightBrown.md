---
source: src/graphics/graphics_init.c
---
**Definition.** `LightBrown` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightBrown", rgb_color(0.94, 0.91, 0.88))` in
`graphics_init` (`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the
literal `RGBColor[0.94, 0.91, 0.88]` and marks it `Protected`. The docstring is kept
centrally in `src/info.c`.

**Representation.** `LightBrown` evaluates to `RGBColor[0.94, 0.91, 0.88]`, so
`FullForm[LightBrown]` is `RGBColor[0.94, 0.91, 0.88]` and `Head[LightBrown]` is
`RGBColor`. Its components are machine reals (none is exactly `0` or `1`). It is a
*separate* constant from `Brown`, not a computed tint of it — the `Light*` family is its
own set of pale pastel entries (`LightRed`, `LightGreen`, `LightBlue`, …).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightBrown, Disk[]}]`, `PlotStyle -> LightBrown`,
`Directive[LightBrown, Thick]` — and as a pale background-style fill. It carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.
