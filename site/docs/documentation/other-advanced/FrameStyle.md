# FrameStyle

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`FrameStyle`**

is an option for Graphics and Plot that specifies the style of the frame box, its ticks and its labels.

<details>
<summary>Notes</summary>

FrameStyle -\> RGBColor\[...\] or GrayLevel\[...\] sets the colour. The default is a neutral gray.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare option symbol, not a function

```mathematica
In[1]:= Head[FrameStyle]
Out[1]= Symbol
```

The named colour resolves to GrayLevel[0]

```mathematica
In[2]:= FullForm[FrameStyle -> Black]
Out[2]= Rule[FrameStyle, GrayLevel[0]]
```

Accepted, still a Graphics

```mathematica
In[3]:= Head[Graphics[{Line[{{0, 0}, {1, 1}}]}, FrameStyle -> Red]]
Out[3]= Graphics
```

## Implementation notes

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

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`FrameStyle` is an option for `Graphics` and `Plot` that styles the frame box, its ticks
and its labels. It is an inert option-name symbol — no builtin, no DownValues, not even
`Protected` — so it does nothing on its own; the frame-drawing code reads it from a call's
option list.

Its value is a style directive, usually a colour. Because named colours are OwnValues,
`FrameStyle -> Black` is stored as `Rule[FrameStyle, GrayLevel[0]]` and `FrameStyle -> Red`
as `Rule[FrameStyle, RGBColor[1, 0, 0]]`. The option has effect only when a `Frame` is
actually drawn, and never alters the plotted data. The default frame colour is a neutral
gray.
