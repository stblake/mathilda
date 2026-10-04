# ContourLabels

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ContourLabels`**

ContourPlot option: True draws z-value text labels at the first visible point of each contour level. Default False.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

An inert, Protected option keyword

```mathematica
In[1]:= Attributes[ContourLabels]
Out[1]= {Protected}
```

Label each level

```mathematica
In[2]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, ContourLabels -> True]]
Out[2]= Graphics
```

## Implementation notes

**Definition.** `ContourLabels` is an **option symbol** for `ContourPlot` that controls
whether each contour level is annotated with its `z`-value. It has no builtin and no
value — an inert, `Protected` keyword read out of `ContourPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_ContourLabels`). `ContourPlot`'s
option handler (`src/graphics/contourplot.c`) reads a `ContourLabels -> b` rule: `True`
draws `z`-value text at the first visible point of each contour level, and the default
`False` draws no labels.

**Usage & limits.** `Protected`. Meaningful only inside a `ContourPlot` call; on its own
it evaluates to itself, and it is read directly from the option sequence rather than
from an `Options[ContourPlot]` default. Related keywords are `Contours`, `ContourStyle`
and `ContourShading`.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/contourplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/contourplot.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ContourLabels` is a `ContourPlot` option: `True` writes each contour's `z`-value as text
at the first visible point of that level, and the default `False` draws no labels.

It is an inert, `Protected` keyword with no value of its own — read straight from the
`ContourPlot` option sequence rather than from a registered `Options` default, so on its
own it just evaluates to itself. Related keywords are `Contours`, `ContourStyle` and
`ContourShading`.
