---
source: src/info.c
---
**Definition.** `ImageSize` is an **option name** for `Graphics` and `Plot` giving the
overall size of the image to display. It is a bare option symbol, not a function: no
builtin, no DownValues, and not even the `Protected` attribute — its only registration is
the docstring set in `info_init` (`src/info.c`). Its meaning lives in the graphics back end
(`src/graphics/`), which reads it from a call's option list when sizing the window or
export.

**Representation.** `ImageSize` stays an inert `EXPR_SYMBOL` and appears only on the left
of a rule, so `FullForm[ImageSize -> 400]` is `Rule[ImageSize, 400]`. The value is either a
width `w` in pixels (the height then following from `AspectRatio`) or a `{w, h}` pair
fixing both. The default width is 800.

**Usage & limits.** Supplied as `Graphics[prims, ImageSize -> 400]` or
`Plot[f, {x, a, b}, ImageSize -> {600, 400}]`, it is consumed only at render/export time
and never changes the plotted data — a call carrying it still reduces to an ordinary
`Graphics[...]` object. It is not among the rules returned by `Options[Plot]`, since `Plot`
forwards it to the renderer rather than storing it as a plotting default.
