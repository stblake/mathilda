# $PlotLegendData

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$PlotLegendData[{color1, label1}, ...]`**

Internal Plot metadata read by the renderer to draw a legend box. Not intended for direct use.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Plot embeds one legend-data node

```mathematica
In[1]:= Count[Plot[Sin[x], {x, 0, Pi}, PlotLegends -> {"sin"}], _$PlotLegendData, Infinity]
Out[1]= 1
```

Its {colour, label} entries

```mathematica
In[2]:= FullForm[$PlotLegendData[{RGBColor[1, 0, 0], "sin"}, {RGBColor[0, 0, 1], "cos"}]]
Out[2]= $PlotLegendData[List[RGBColor[1, 0, 0], "sin"], List[RGBColor[0, 0, 1], "cos"]]
```

An inert metadata head

```mathematica
In[3]:= MatchQ[$PlotLegendData[{RGBColor[1, 0, 0], "sin"}], _$PlotLegendData]
Out[3]= True
```

## Implementation notes

**Definition.** `$PlotLegendData` is an **internal metadata head** — not a user-facing
function. It carries the legend of a plot inside the option list of a `Graphics[...]`
object: `$PlotLegendData[{colour, label}, {colour, label}, ...]`, one entry per curve.
It is built by `Plot` (and friends) when `PlotLegends` is given, and read back by the
renderer when it draws the legend box.

**Representation.** A bare `EXPR_SYMBOL` head (interned `SYM_PlotLegendData`) with the
attributes `{HoldAll, Protected}`. `plot.c`'s `build_legend_meta` constructs it from the
`PlotLegends` option and appends it to the `Graphics[...]` option sequence;
`render.c`'s `find_legend_data` locates it there, and `Show[]` concatenates the
`$PlotLegendData` entries of several graphics into one combined legend (while dropping
`$PlotResample`). `HoldAll` keeps the stored `{colour, label}` pairs from re-evaluating.

**Usage & limits.** Internal plumbing: a user never writes `$PlotLegendData` directly,
and it is meaningful only as an item in a graphics option list. Evaluated on its own it
stays symbolic (it has no evaluation rule). It is documented here for completeness
because it is reachable from the symbol table (`Names["*"]`).

*Internal symbol — shown below only structurally; it has no standalone mathematical
example.*

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/graphics/render.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/render.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`$PlotLegendData[{colour, label}, ...]` is an **internal metadata head** — not something
you write by hand. When `Plot` (or a related head) is given `PlotLegends`, it attaches a
`$PlotLegendData` node to the `Graphics[...]` option list, one `{colour, label}` entry
per curve; the renderer reads it back to draw the legend box, and `Show[]` concatenates
the nodes of several graphics into one combined legend.

It carries `{HoldAll, Protected}`, so its stored pairs do not re-evaluate. On its own it
has no evaluation rule and stays symbolic; it is listed here only because it is reachable
from the symbol table. The examples above are structural — this symbol has no standalone
mathematical meaning.
