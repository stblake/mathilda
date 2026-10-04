### Worked examples

```mathematica
In[1]:= Attributes[Frame]  (* a bare option keyword -- no attributes, not even Protected *)
```

```mathematica
In[1]:= Head[Plot[Sin[x], {x, 0, Pi}, Frame -> True]]  (* a framed plot is still a Graphics object *)
```

### Notes

`Frame` is an option for `Graphics` and the plotting heads that decides whether a boxed
frame with ticks and labels is drawn around the plot. `Frame -> True` boxes all four
edges, `Frame -> False` (or `None`) draws none, and `{{left, right}, {bottom, top}}`
toggles each edge individually. In `Plot` a frame takes the place of the default `Axes`;
some heads (such as `ArrayPlot` and `ComplexPlot`) default to `Frame -> True`.

Unlike the other graphics option keywords, `Frame` is **not** `Protected` — its
`Attributes` are empty — so it can in principle be reassigned. On its own it evaluates to
itself; it is meaningful only inside a graphics or plot call. Companion keywords are
`FrameTicks`, `FrameLabel`, `Axes` and `AxesLabel`.
