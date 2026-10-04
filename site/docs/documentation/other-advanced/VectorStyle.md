# VectorStyle

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`VectorStyle`**

VectorPlot option: style directive(s) (RGBColor, Thickness, …) applied globally to all arrows. Overrides per-arrow ColorFunction.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare VectorPlot option symbol

```mathematica
In[1]:= Head[VectorStyle]
Out[1]= Symbol
```

Inert and Protected -- no builtin behaviour

```mathematica
In[2]:= Attributes[VectorStyle]
Out[2]= {Protected}
```

The named colour resolves to RGBColor[1, 0, 0]

```mathematica
In[3]:= FullForm[VectorStyle -> Red]
Out[3]= Rule[VectorStyle, RGBColor[1, 0, 0]]
```

## Implementation notes

**Definition.** `VectorStyle` is an **option name** for `VectorPlot` giving the style
directive(s) applied globally to every arrow. It is a bare option symbol, not a function:
`graphics_init` (`src/graphics/graphics_init.c`) stamps it `Protected` and sets its
docstring, but it has no builtin and no DownValues. The `VectorPlot` renderer reads it from
the call's option list.

**Representation.** `VectorStyle` stays an inert, `Protected` `EXPR_SYMBOL`, appearing only
on the left of a rule; its value is one or more style directives (`RGBColor`, `Thickness`,
…). Because named colours are OwnValues, `FullForm[VectorStyle -> Red]` resolves to
`Rule[VectorStyle, RGBColor[1, 0, 0]]`.

**Usage & limits.** Supplied as
`VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorStyle -> Red]`, it applies the given
directive to all arrows and overrides any per-arrow `ColorFunction`. It is consumed at
render time and affects only the arrows' appearance, never the sampled field. It is
specific to `VectorPlot`, not a general `Graphics` directive.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`VectorStyle` is an option for `VectorPlot` that applies style directive(s) globally to
every arrow. It is an inert, `Protected` option-name symbol — no builtin, no DownValues —
so it does nothing on its own; the `VectorPlot` renderer reads it from the call's option
list.

Its value is one or more directives (`RGBColor`, `Thickness`, …); a named colour resolves,
so `VectorStyle -> Red` is stored as `Rule[VectorStyle, RGBColor[1, 0, 0]]`. It overrides
any per-arrow `ColorFunction` and affects only the arrows' appearance, never the sampled
field.
