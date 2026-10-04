---
source: src/graphics/graphics_init.c
---
**Definition.** `Black` is one of Mathilda's named colour constants. It is bound as an
**OwnValue** — `register_color("Black", gray_color(0))` in `graphics_init`
(`src/graphics/graphics_init.c`) — so the bare symbol `Black` evaluates to the literal
`GrayLevel[0]` and is marked `Protected`. Unlike the saturated primaries it uses
`GrayLevel`, not `RGBColor`, matching Mathematica's own `InputForm[Black]`. The docstring
is kept centrally in `src/info.c`.

**Representation.** `Black` evaluates to `GrayLevel[0]`, so `FullForm[Black]` is
`GrayLevel[0]` and `Head[Black]` is `GrayLevel`. The single component is an exact integer
`0` (the `color_component` helper prints an exact `0`/`1` as an integer). `Black`, `White`
(`GrayLevel[1]`), `Gray` (`GrayLevel[0.5]`) and `LightGray` (`GrayLevel[0.85]`) are the
greyscale family.

**Usage & limits.** As a resolved directive it is accepted wherever a colour is —
`Graphics[{Black, Line[...]}]`, `PlotStyle -> Black`, `FrameStyle -> Black`,
`Directive[Black, Thick]`. A renderer only ever has to understand the long `GrayLevel[...]`
form. It carries no opacity; use `Opacity[a]` or `GrayLevel[g, a]` for a translucent black.
