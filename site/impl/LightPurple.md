---
source: src/graphics/graphics_init.c
---
**Definition.** `LightPurple` is a named colour constant standing for `RGBColor[0.94, 0.88, 0.94]` — a pale purple tint
in the additive (screen) RGB model. It is **not** an inert head: `register_color`
(`graphics_init.c`) installs it as an **OwnValue**, so `LightPurple` *evaluates* to its
`RGBColor[0.94, 0.88, 0.94]` literal, and the `Protected` attribute is set on the same symbol. The
docstring is kept centrally in `info.c`, mirroring `CMYKColor` and the other named
colours.

**Representation.** The OwnValue is built by `rgb_color(...)`, which runs each
channel through `color_component`: a component exactly `0` or `1` is emitted as a
plain `Integer` (matching Mathematica's `FullForm`), any other value as a `Real`.
Because `LightPurple` reduces to a real colour literal before it ever reaches a renderer,
it resolves anywhere `RGBColor`/`GrayLevel`/`Hue` do — a primitive-list directive,
a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction` return.

**Usage & limits.** A style directive only: it sets the colour of the graphics
primitives that follow it and carries no numeric meaning elsewhere — `D`, `N` and
arithmetic see the `RGBColor[0.94, 0.88, 0.94]` triple, not a scalar. The named colours are defined
together as one table in `graphics_init.c`; `LightPurple` is a single row of it.
