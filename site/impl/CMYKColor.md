---
source: src/graphics/render.c
---
**Algorithm.** `CMYKColor` is an inert, `Protected` **style directive** for the
subtractive (print) colour model — it has no builtin and no OwnValue (only the
attribute bit is set in `graphics_init.c`; its docstring lives centrally in
`info.c`), so it stays symbolic and is interpreted only at render time, exactly
like `RGBColor`/`GrayLevel`/`Hue`. The forms are `CMYKColor[c, m, y, k]`,
`CMYKColor[c, m, y]` (`k = 0`), `CMYKColor[c, m, y, k, a]` (with opacity), and the
list forms `CMYKColor[{c, m, y, k}]` / `CMYKColor[{c, m, y, k, a}]`. At render
time `rgba_from_cmyk`/`cmyk_to_rgb` (`render.c`) clip every component to `[0,1]`
(`clip01`) and convert with `w = 1 - k`, `r = (1 - c) w`, `g = (1 - m) w`,
`b = (1 - y) w`. It is recognised as a colour directive by the interned-pointer
test `h == SYM_CMYKColor` in `resolve_color` wherever `RGBColor`/`GrayLevel`/`Hue`
are — as a primitive-list directive, a `PlotStyle`/`Background`/`FrameStyle`
value, or a `ColorFunction` return value. The **same** conversion is duplicated in
`graphics_json.c` for the JSON/notebook-frontend path, so the CLI Raylib renderer
and the notebook agree.

**Data structures.** No dedicated CMYK struct: it is converted straight to an
`RGBA8 {unsigned char r, g, b, a}` at parse/render time (and thence to Raylib's
`Color` via `to_raylib`), or serialised to an `rgba(...)` string in the JSON
path.

**Complexity / limits.** `O(1)` per colour. Components and opacity outside
`[0,1]` are clipped rather than rejected. Because it is inert, `CMYKColor[...]`
never evaluates to an `RGBColor` at the language level — the conversion happens
only inside the two renderers.
