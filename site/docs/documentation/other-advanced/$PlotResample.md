# $PlotResample

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`$PlotResample[var, {f...}, {plotPoints, maxRecursion, maxPlotPoints, mesh, regionFunction, exclusions, colorFunction, colorFunctionScaling, filling, fillingStyle}]`**

Internal Plot metadata used by the renderer to re-sample curves at the current zoom. Not intended for direct use.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare internal metadata symbol

```mathematica
In[1]:= Head[$PlotResample]
Out[1]= Symbol
```

HoldAll keeps the plotted bodies unevaluated

```mathematica
In[2]:= Attributes[$PlotResample]
Out[2]= {HoldAll, Protected}
```

No standalone meaning -- stays unevaluated

```mathematica
In[3]:= $PlotResample[x, {Sin[x]}, {}]
Out[3]= $PlotResample[x, {Sin[x]}, {}]
```

## Implementation notes

**Definition.** `$PlotResample` is an **internal metadata head** that `Plot` embeds inside
the `Graphics` object it returns, so the renderer can re-sample the curves when the view is
zoomed. It is not a user function: `graphics_init` (`src/graphics/graphics_init.c`) stamps
it `HoldAll | Protected` and sets a docstring that says outright it is not intended for
direct use; the node is built by `make_resample_meta` in `src/graphics/plot.c` and consumed
by the renderer in `src/graphics/show.c`.

**Representation.** `$PlotResample` is a bare `EXPR_SYMBOL` (`Head[$PlotResample]` is
`Symbol`, `Attributes[$PlotResample]` is `{HoldAll, Protected}`). In use it heads a
three-argument form `$PlotResample[var, {bodies}, {opts...}]`, where `bodies` are the
unevaluated plotted expressions and `opts` is the packed option vector (`plotPoints`,
`maxRecursion`, `maxPlotPoints`, mesh, region function, exclusions, colour function and its
scaling, filling, filling style). `HoldAll` is what keeps the bodies and the held option
values unevaluated through the surrounding `Graphics` re-evaluation.

**Usage & limits.** It exists only as a carrier: `show.c` recognises the head while drawing
and, on a zoom, re-samples the stored bodies over the new coordinate window instead of
rescaling a fixed polyline. It has no standalone value — written by hand it simply stays
unevaluated — and is an implementation detail of `Plot`'s zoomable output.

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/graphics/plot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/plot.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`$PlotResample` is internal `Plot` metadata, not a user-facing function. `Plot` embeds a
`$PlotResample[var, {bodies}, {opts...}]` node inside the `Graphics` object it returns so
the renderer can re-sample the curves when the view is zoomed, rather than rescaling a fixed
polyline.

Its `HoldAll` attribute keeps the plotted bodies and held option values unevaluated through
the surrounding `Graphics` re-evaluation. It carries no standalone value — written by hand
it simply stays unevaluated, as above — and should be treated as an implementation detail of
`Plot`'s zoomable output.
