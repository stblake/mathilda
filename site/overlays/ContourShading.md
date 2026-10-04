### Worked examples

```mathematica
In[1]:= Head[ContourShading]  (* a bare ContourPlot option symbol *)
```

```mathematica
In[1]:= Attributes[ContourShading]  (* inert and Protected -- no builtin behaviour *)
```

```mathematica
In[1]:= FullForm[ContourShading -> False]  (* it only sits on the left of a rule *)
```

### Notes

`ContourShading` is an option for `ContourPlot` that controls whether the cells between
contour lines are filled. It is an inert, `Protected` option-name symbol — no builtin, no
DownValues — so it does nothing on its own; the `ContourPlot` sampler reads it from the
call's option list.

`ContourShading -> True` fills each cell by its `z` value, using the supplied
`ColorFunction` or the built-in blue-cyan-yellow-red thermal ramp; `False` draws contour
lines only; `Automatic` shades only when a `ColorFunction` is given. It affects the fill
only, never the sampled data, and is specific to `ContourPlot`.
