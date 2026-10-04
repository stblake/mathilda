### Worked examples

```mathematica
In[1]:= Head[AspectRatio]  (* a bare option symbol, not a function *)
```

```mathematica
In[1]:= FullForm[AspectRatio -> 1]  (* it only ever sits on the left of a rule *)
```

```mathematica
In[1]:= Options[Plot, AspectRatio]  (* Plot's stored default is Automatic *)
```

```mathematica
In[1]:= Head[Graphics[{Disk[]}, AspectRatio -> 1, ImageSize -> 300]]  (* accepted, still a Graphics *)
```

### Notes

`AspectRatio` is an option for `Graphics` and `Plot` that fixes the ratio of rendered
height to width. It is an inert option-name symbol — no builtin, no DownValues, not even
`Protected` — so it does nothing on its own; the graphics back end reads it out of a call's
option list.

`AspectRatio -> r` sets an explicit height-to-width ratio; `AspectRatio -> Automatic` takes
the true geometry from the coordinate values; `AspectRatio -> Full` stretches to fill the
region. `Plot` reports a stored default of `Automatic` and otherwise lays a standalone
curve out at `1/GoldenRatio`. The option changes the displayed shape only, never the
plotted data, so a call carrying it still evaluates to an ordinary `Graphics[...]` object.
