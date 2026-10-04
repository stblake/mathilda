### Worked examples

```mathematica
In[1]:= Count[Plot[Sin[x], {x, 0, Pi}, PlotLegends -> {"sin"}], _$PlotLegendData, Infinity]  (* Plot embeds one legend-data node *)
```

```mathematica
In[1]:= FullForm[$PlotLegendData[{RGBColor[1, 0, 0], "sin"}, {RGBColor[0, 0, 1], "cos"}]]  (* its {colour, label} entries *)
```

```mathematica
In[1]:= MatchQ[$PlotLegendData[{RGBColor[1, 0, 0], "sin"}], _$PlotLegendData]  (* an inert metadata head *)
```

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
