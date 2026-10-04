# ContourShading

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ContourShading`**

ContourPlot option: True/False/Automatic — fills cells by z value using ColorFunction or the built-in blue-cyan-yellow-red thermal ramp. Automatic enables shading only when ColorFunction is set.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare ContourPlot option symbol

```mathematica
In[1]:= Head[ContourShading]
Out[1]= Symbol
```

Inert and Protected -- no builtin behaviour

```mathematica
In[2]:= Attributes[ContourShading]
Out[2]= {Protected}
```

It only sits on the left of a rule

```mathematica
In[3]:= FullForm[ContourShading -> False]
Out[3]= Rule[ContourShading, False]
```

## Implementation notes

**Definition.** `ContourShading` is an **option name** for `ContourPlot` that controls
whether the grid cells between contour lines are filled. It is a bare option symbol, not a
function: `graphics_init` (`src/graphics/graphics_init.c`) stamps it `Protected` and sets
its docstring, but it has no builtin and no DownValues. The `ContourPlot` sampler reads it
from the call's option list.

**Representation.** `ContourShading` stays an inert, `Protected` `EXPR_SYMBOL`, appearing
only on the left of a rule; `FullForm[ContourShading -> False]` is
`Rule[ContourShading, False]`. The value is `True`, `False`, or `Automatic`.

**Usage & limits.** With `ContourShading -> True` each cell is filled by its `z` value,
using the supplied `ColorFunction` or the built-in blue-cyan-yellow-red thermal ramp;
`False` draws only the contour lines; `Automatic` enables shading only when a
`ColorFunction` is set. It is consumed at render time and affects only the fill, never the
sampled data. It is specific to `ContourPlot`, not a general `Graphics` directive.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ContourShading` is an option for `ContourPlot` that controls whether the cells between
contour lines are filled. It is an inert, `Protected` option-name symbol — no builtin, no
DownValues — so it does nothing on its own; the `ContourPlot` sampler reads it from the
call's option list.

`ContourShading -> True` fills each cell by its `z` value, using the supplied
`ColorFunction` or the built-in blue-cyan-yellow-red thermal ramp; `False` draws contour
lines only; `Automatic` shades only when a `ColorFunction` is given. It affects the fill
only, never the sampled data, and is specific to `ContourPlot`.
