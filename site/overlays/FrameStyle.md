### Worked examples

```mathematica
In[1]:= Head[FrameStyle]  (* a bare option symbol, not a function *)
```

```mathematica
In[1]:= FullForm[FrameStyle -> Black]  (* the named colour resolves to GrayLevel[0] *)
```

```mathematica
In[1]:= Head[Graphics[{Line[{{0, 0}, {1, 1}}]}, FrameStyle -> Red]]  (* accepted, still a Graphics *)
```

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
