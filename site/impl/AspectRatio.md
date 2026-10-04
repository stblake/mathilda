---
source: src/info.c
---
**Definition.** `AspectRatio` is an **option name** for `Graphics` and `Plot` giving the
ratio of rendered height to width. It is a bare option symbol, not a function: it has no
builtin, no DownValues and (unusually) not even the `Protected` attribute — its only
registration is the docstring set in `info_init` (`src/info.c`). Its meaning lives entirely
in the graphics back end (`src/graphics/`), which reads the option out of the option list
a `Graphics`/`Plot` call carries.

**Representation.** `AspectRatio` stays an inert `EXPR_SYMBOL`; it appears only on the
left of a rule, so `FullForm[AspectRatio -> 1]` is `Rule[AspectRatio, 1]`. The value may
be an explicit height-to-width real, `Automatic` (derive the ratio from the actual
coordinate values — true geometry), or `Full` (stretch to fill the enclosing region).
`Options[Plot, AspectRatio]` reports `Plot`'s stored default, `Automatic`; a standalone
`Plot` curve otherwise lays out at `1/GoldenRatio`.

**Usage & limits.** Supplied as `Graphics[prims, AspectRatio -> r]` or
`Plot[f, {x, a, b}, AspectRatio -> r]`, it is consumed by the renderer when the window or
export is laid out; it changes the displayed shape, never the plotted data, so a call that
accepts it still reduces to an ordinary `Graphics[...]` object. Because the symbol is not
`Protected`, it could in principle be assigned a value — it is intended only as an option
tag.
