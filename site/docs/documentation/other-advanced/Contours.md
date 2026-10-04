# Contours

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Contours`**

ContourPlot option: integer count of auto-levels or explicit list of contour values.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

An inert, Protected option keyword

```mathematica
In[1]:= Attributes[Contours]
Out[1]= {Protected}
```

5 auto-levels

```mathematica
In[2]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, Contours -> 5]]
Out[2]= Graphics
```

## Implementation notes

**Definition.** `Contours` is an **option symbol** for `ContourPlot` that selects the
contour levels. It has no builtin and no value — it is an inert, `Protected` keyword
read out of `ContourPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Contours`). `ContourPlot`'s
option handler (`src/graphics/contourplot.c`) reads a `Contours -> spec` rule where
`spec` is either an integer `n` (draw `n` automatically chosen, evenly spaced levels)
or an explicit list of `z`-values at which to draw contours. The chosen levels drive
the marching-squares sampler and, when shading is on, the colour ramp.

**Usage & limits.** `Protected`. Meaningful only inside a `ContourPlot` call; on its own
it evaluates to itself. It is read directly from the option sequence rather than from a
registered default in `Options[ContourPlot]`. Related keywords are `ContourStyle`,
`ContourLabels` and `ContourShading`.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/contourplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/contourplot.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Contours` is a `ContourPlot` option selecting the contour levels: an integer `n` draws
`n` automatically chosen, evenly spaced levels, and an explicit list draws contours at
exactly those `z`-values.

It is an inert, `Protected` keyword with no value of its own — read straight from the
`ContourPlot` option sequence rather than from a registered `Options` default, so on its
own it just evaluates to itself. Related keywords are `ContourStyle`, `ContourLabels`
and `ContourShading`.
