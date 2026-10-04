### Worked examples

```mathematica
In[1]:= Count[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, PlotLegends -> Automatic], _$StreamColorBar, Infinity]  (* ContourPlot embeds one colour-bar node *)
```

```mathematica
In[1]:= Head[$StreamColorBar[-Pi, Pi, "Cyclic"]]  (* an inert metadata head *)
```

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
