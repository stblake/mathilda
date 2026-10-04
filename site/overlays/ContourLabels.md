### Worked examples

```mathematica
In[1]:= Attributes[ContourLabels]  (* an inert, Protected option keyword *)
```

```mathematica
In[1]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, ContourLabels -> True]]  (* label each level *)
```

### Notes

`ContourLabels` is a `ContourPlot` option: `True` writes each contour's `z`-value as text
at the first visible point of that level, and the default `False` draws no labels.

It is an inert, `Protected` keyword with no value of its own — read straight from the
`ContourPlot` option sequence rather than from a registered `Options` default, so on its
own it just evaluates to itself. Related keywords are `Contours`, `ContourStyle` and
`ContourShading`.
