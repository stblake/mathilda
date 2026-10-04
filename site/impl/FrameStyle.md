---
source: src/info.c
---
**Definition.** `FrameStyle` is an **option name** for `Graphics` and `Plot` giving the
style of the frame box, its ticks and its labels. It is a bare option symbol, not a
function: no builtin, no DownValues, not even the `Protected` attribute — its only
registration is the docstring in `info_init` (`src/info.c`). The frame-drawing code in the
graphics back end (`src/graphics/`) reads it from a call's option list, alongside the
related `Frame` and `FrameTicks` options.

**Representation.** `FrameStyle` stays an inert `EXPR_SYMBOL`, appearing only on the left
of a rule; its value is a style directive — typically a colour such as `RGBColor[...]` or
`GrayLevel[...]`. Because a named colour is an OwnValue, `FullForm[FrameStyle -> Black]`
resolves to `Rule[FrameStyle, GrayLevel[0]]` and `FrameStyle -> Red` to
`Rule[FrameStyle, RGBColor[1, 0, 0]]`. The default frame colour is a neutral gray.

**Usage & limits.** Supplied as `Plot[f, {x, a, b}, Frame -> True, FrameStyle -> Red]` or
inside a `Graphics[...]`, it is consumed only when the frame is drawn and has no effect
unless a frame is present. It never changes the plotted data, so a call carrying it still
reduces to an ordinary `Graphics[...]` object.
