# VectorScale

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`VectorScale`**

VectorPlot option: controls arrow length. Automatic (default): all arrows drawn at equal length (direction only) None: length proportional to field magnitude real f: arrow length = f × grid spacing

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert, Protected option keyword

```mathematica
In[1]:= Attributes[VectorScale]
Out[1]= {Protected}
```

No value of its own: evaluates to itself

```mathematica
In[2]:= VectorScale
Out[2]= VectorScale
```

Read inside a VectorPlot

```mathematica
In[3]:= Head[VectorPlot[{1, x}, {x, -1, 1}, {y, -1, 1}, VectorScale -> 0.5]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `VectorScale` is an **option symbol** for `VectorPlot` that controls how
arrow length encodes field magnitude. It has no builtin and no value of its own — it is
an inert keyword (`Protected`) read out of `VectorPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_VectorScale`). `VectorPlot`'s
option parser (`src/graphics/vectorplot.c`) evaluates the right-hand side of a
`VectorScale -> v` rule and sets an internal `scale_mode`: `Automatic` (the default) →
`VS_AUTOMATIC`, all arrows drawn at equal display length so only direction is shown;
`None` → `VS_NONE`, length proportional to field magnitude; a positive real `f` →
`VS_FRACTION` with arrow length `= f × grid spacing`. Any other value falls back to
`Automatic`.

**Usage & limits.** `Protected`. It carries meaning only inside a `VectorPlot` call; on
its own it just evaluates to itself. It has no entry in `Options[VectorPlot]`'s default
list — it is read directly from the given option sequence — so `SetOptions` on it has no
effect. Companion `VectorPlot` options are `VectorPoints` (seed-grid density) and
`VectorStyle` (global style directives).

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/vectorplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/vectorplot.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`VectorScale` is a `VectorPlot` option controlling how arrow length encodes the field.
`Automatic` (the default) draws every arrow at equal length, showing direction only;
`None` makes length proportional to the field magnitude; a positive real `f` sets arrow
length to `f ×` the grid spacing.

It is an inert, `Protected` keyword with no value of its own, read straight from the
`VectorPlot` option sequence (it is not in `Options[VectorPlot]`'s default list, so
`SetOptions` on it has no effect). The companion keywords are `VectorPoints` (seed-grid
density) and `VectorStyle` (global style directives).
