### Worked examples

```mathematica
In[1]:= Head[ImageSize]  (* a bare option symbol, not a function *)
```

```mathematica
In[1]:= FullForm[ImageSize -> 400]  (* it only ever sits on the left of a rule *)
```

```mathematica
In[1]:= Head[Graphics[{Disk[]}, ImageSize -> {600, 400}]]  (* accepted, still a Graphics *)
```

### Notes

`ImageSize` is an option for `Graphics` and `Plot` that sets the overall display size. It
is an inert option-name symbol — no builtin, no DownValues, not even `Protected` — so it
does nothing on its own; the graphics back end reads it from a call's option list when it
sizes the window or an export.

`ImageSize -> w` sets the width in pixels, with the height following from `AspectRatio`;
`ImageSize -> {w, h}` fixes both. The default width is 800. The option affects only the
rendered size, never the plotted data, so a call carrying it still evaluates to an ordinary
`Graphics[...]` object.
