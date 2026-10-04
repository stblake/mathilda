---
source: src/graphics/graphics_init.c
---
**Definition.** `LightGray` is one of Mathilda's **named colour constants** — a pale gray. Like
Mathematica's own named colours it is not an inert head but a symbol whose OwnValue is
a `GrayLevel` literal, so it resolves to a real colour directive anywhere one is expected:
`LightGray` evaluates to `GrayLevel[0.85]`.

**Representation.** `register_color("LightGray", ...)` in `src/graphics/graphics_init.c`
installs the OwnValue `LightGray -> GrayLevel[0.85]` and marks the symbol `Protected`; the
component builder prints an exact `0` or `1` as a plain integer and anything else as
a real, matching Mathematica's own `InputForm`. The docstring lives centrally in
`src/info.c`. Because it is an OwnValue (not an inert head), a non-held graphics
argument or an evaluated-once plot option sees the expanded `GrayLevel[0.85]`, exactly as a
literal colour would.

**Usage & limits.** `Protected`. It is a value, not a function: feeding it a wrong
"argument" is meaningless. Use it as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`. The grey-scale constants
(`Black`, `White`, `Gray`, `LightGray`) expand to `GrayLevel`; the rest to
`RGBColor`.
