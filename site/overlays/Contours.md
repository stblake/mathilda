### Worked examples

```mathematica
In[1]:= Attributes[Contours]  (* an inert, Protected option keyword *)
```

```mathematica
In[1]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, Contours -> 5]]  (* 5 auto-levels *)
```

### Notes

`Contours` is a `ContourPlot` option selecting the contour levels: an integer `n` draws
`n` automatically chosen, evenly spaced levels, and an explicit list draws contours at
exactly those `z`-values.

It is an inert, `Protected` keyword with no value of its own — read straight from the
`ContourPlot` option sequence rather than from a registered `Options` default, so on its
own it just evaluates to itself. Related keywords are `ContourStyle`, `ContourLabels`
and `ContourShading`.
