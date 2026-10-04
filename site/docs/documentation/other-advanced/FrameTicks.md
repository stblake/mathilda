# FrameTicks

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`FrameTicks`**

is an option for Graphics and Plot that specifies the tick marks on the edges of a frame.

<details>
<summary>Notes</summary>

FrameTicks -\> Automatic (the default) draws major and minor ticks with labels on every drawn frame edge; FrameTicks -\> None keeps the frame box but draws no ticks; the {{left, right}, {bottom, top}} form selects Automatic or None per edge.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= FrameTicks -> None
Out[1]= FrameTicks -> None
```

Per-edge forms pass through unevaluated too

```mathematica
In[2]:= FrameTicks -> {Automatic, None}
Out[2]= FrameTicks -> {Automatic, None}
```

Accepted by Plot, which builds a Graphics

```mathematica
In[3]:= Head[Plot[Sin[x], {x, 0, 3}, FrameTicks -> None]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `FrameTicks` is an option for `Graphics` and `Plot` that specifies
the tick marks on the edges of a frame. `FrameTicks -> Automatic` (the default)
draws major and minor ticks with labels on every drawn frame edge;
`FrameTicks -> None` keeps the frame box but draws no ticks; the
`{{left, right}, {bottom, top}}` form selects `Automatic` or `None` per edge. It is
an inert option keyword — no builtin and no value of its own (it carries no
attributes), declared with its docstring in `info.c`.

**Representation.** A bare `EXPR_SYMBOL`. It is read as the left of a
`FrameTicks -> value` rule by the 2D renderer in `src/graphics/render.c`: the option
loop matches `name == SYM_FrameTicks` and interprets the value — `Automatic`/`True`
turns ticks on for every drawn edge, `None`/`False` keeps the frame box but
suppresses ticks and labels, and a `{{l, r}, {b, t}}` nesting is read per edge.

**Usage & limits.** Meaningful only inside a `Graphics`/`Plot` call that draws a
frame; with `Frame -> False` there are no frame edges to tick. The symbol itself
performs no computation — all behaviour is in the renderer's reading of the rule.

**Attributes:** none registered.

## References

- Source: [`src/graphics/render.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/render.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`FrameTicks` is an option for `Graphics` and `Plot` controlling the tick marks on a
frame's edges: `Automatic` (the default) draws labelled major and minor ticks on
every drawn edge, `None` keeps the frame box but no ticks, and a
`{{left, right}, {bottom, top}}` form selects `Automatic` or `None` per edge. It is
an inert keyword with no value of its own; all behaviour lives in the 2D renderer's
reading of the `FrameTicks -> value` rule. It only bites when a frame is actually
drawn.
