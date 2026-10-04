### Worked examples

```mathematica
In[1]:= Head[$PlotResample]  (* a bare internal metadata symbol *)
```

```mathematica
In[1]:= Attributes[$PlotResample]  (* HoldAll keeps the plotted bodies unevaluated *)
```

```mathematica
In[1]:= $PlotResample[x, {Sin[x]}, {}]  (* no standalone meaning -- stays unevaluated *)
```

### Notes

`$PlotResample` is internal `Plot` metadata, not a user-facing function. `Plot` embeds a
`$PlotResample[var, {bodies}, {opts...}]` node inside the `Graphics` object it returns so
the renderer can re-sample the curves when the view is zoomed, rather than rescaling a fixed
polyline.

Its `HoldAll` attribute keeps the plotted bodies and held option values unevaluated through
the surrounding `Graphics` re-evaluation. It carries no standalone value — written by hand
it simply stays unevaluated, as above — and should be treated as an implementation detail of
`Plot`'s zoomable output.
