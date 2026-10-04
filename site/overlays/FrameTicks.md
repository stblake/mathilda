### Worked examples

```mathematica
In[1]:= FrameTicks -> None  (* an inert option keyword: it just names the left of a rule *)
```

```mathematica
In[1]:= FrameTicks -> {Automatic, None}  (* per-edge forms pass through unevaluated too *)
```

```mathematica
In[1]:= Head[Plot[Sin[x], {x, 0, 3}, FrameTicks -> None]]  (* accepted by Plot, which builds a Graphics *)
```

### Notes

`FrameTicks` is an option for `Graphics` and `Plot` controlling the tick marks on a
frame's edges: `Automatic` (the default) draws labelled major and minor ticks on
every drawn edge, `None` keeps the frame box but no ticks, and a
`{{left, right}, {bottom, top}}` form selects `Automatic` or `None` per edge. It is
an inert keyword with no value of its own; all behaviour lives in the 2D renderer's
reading of the `FrameTicks -> value` rule. It only bites when a frame is actually
drawn.
