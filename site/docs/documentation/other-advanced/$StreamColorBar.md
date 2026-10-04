# $StreamColorBar

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`$StreamColorBar[speed_min, speed_max]`**

Internal StreamPlot metadata. Instructs the renderer to draw a vertical speed color scale bar (dark blue = slow, yellow = fast). Not intended for direct use.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

ContourPlot embeds one colour-bar node

```mathematica
In[1]:= Count[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, PlotLegends -> Automatic], _$StreamColorBar, Infinity]
Out[1]= 1
```

An inert metadata head

```mathematica
In[2]:= Head[$StreamColorBar[-Pi, Pi, "Cyclic"]]
Out[2]= $StreamColorBar
```

## Implementation notes

**Definition.** `$StreamColorBar` is an **internal metadata head** — not a user-facing
function. It records the data needed to draw a colour-scale bar beside a plot:
`$StreamColorBar[min, max]` for a density/contour scale, or
`$StreamColorBar[min, max, colourFnOrName]` for a complex-plot phase wheel. It is
emitted by the plot builders and read back by the renderer.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_StreamColorBar`,
`Protected`). `contourplot.c` appends `$StreamColorBar[zmin, zmax, ...]` to the
`Graphics[...]` option list when `PlotLegends -> Automatic`; `complexplot.c` appends
`$StreamColorBar[-Pi, Pi, cfn_or_"Cyclic"]` for its phase bar; `render.c`'s
`find_stream_colorbar` locates the node in the option list and draws the bar.

**Usage & limits.** Internal plumbing: a user never writes `$StreamColorBar` directly,
and it is meaningful only as an item in a graphics option list. Evaluated on its own it
stays symbolic (no evaluation rule). It is documented here only because it is reachable
from the symbol table (`Names["*"]`).

*Internal symbol — shown below only structurally; it has no standalone mathematical
example.*

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/render.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/render.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`$StreamColorBar[min, max]` (or `$StreamColorBar[min, max, colourFnOrName]`) is an
**internal metadata head** — not something you write by hand. It records the data a
colour-scale bar needs: `ContourPlot` / `DensityPlot` append `$StreamColorBar[zmin,
zmax, ...]` when `PlotLegends -> Automatic`, and `ComplexPlot` appends
`$StreamColorBar[-Pi, Pi, cfn]` for its phase wheel; the renderer reads the node back
and draws the bar.

It is `Protected`, has no evaluation rule, and stays symbolic on its own. It is listed
here only because it is reachable from the symbol table. The examples above are
structural — this symbol has no standalone mathematical meaning.
