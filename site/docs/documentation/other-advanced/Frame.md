# Frame

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Frame`**

is an option for Graphics and Plot that specifies whether to draw a frame with ticks and labels around the plot.

<details>
<summary>Notes</summary>

Frame -\> True boxes all four edges; Frame -\> False (or None) draws no frame; Frame -\> {{left, right}, {bottom, top}} toggles each edge with True or False. In Plot a frame takes the place of the default Axes.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

A bare option keyword -- no attributes, not even Protected

```mathematica
In[1]:= Attributes[Frame]
Out[1]= {}
```

A framed plot is still a Graphics object

```mathematica
In[2]:= Head[Plot[Sin[x], {x, 0, Pi}, Frame -> True]]
Out[2]= Graphics
```

## Implementation notes

**Definition.** `Frame` is an **option symbol** for `Graphics` and the plotting heads
(`Plot`, `ListPlot`, `ContourPlot`, `DensityPlot`, `VectorPlot`, `ArrayPlot`, ...) that
decides whether a boxed frame with ticks and labels is drawn around the plot. It has no
builtin and no value of its own — it is a keyword read out of the option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Frame`). The renderer and each
plot's option parser read a `Frame -> spec` rule: `True` boxes all four edges;
`False`/`None` draws no frame; `{{left, right}, {bottom, top}}` toggles each edge
individually. In `Plot` a frame takes the place of the default `Axes`. Several heads
(e.g. `ArrayPlot`, `ComplexPlot`) default `Frame -> True`.

**Usage & limits.** It carries a docstring from `info.c` but, unlike the other graphics
option keywords, is **not** marked `Protected` (`Attributes[Frame]` is `{}`), so it can
in principle be reassigned. It is meaningful only inside a graphics/plot call; on its
own it evaluates to itself. Companion keywords are `FrameTicks`, `FrameLabel`, `Axes`
and `AxesLabel`.

**Attributes:** none registered.

## References

- Source: [`src/graphics/render.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/render.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

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
