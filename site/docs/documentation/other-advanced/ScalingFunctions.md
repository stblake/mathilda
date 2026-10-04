# ScalingFunctions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ScalingFunctions`**

Option for Plot, ListPlot, DensityPlot, ContourPlot, VectorPlot, StreamPlot: applies a coordinate transform to one or both axes. Forms: ScalingFunctions -\> "Log"         both axes: natural log ScalingFunctions -\> "Log10"        both axes: log base 10 ScalingFunctions -\> "Log2"         both axes: log base 2 ScalingFunctions -\> "Reverse"      both axes: mirror (negate) ScalingFunctions -\> {"Log", None}  x-axis log, y-axis linear ScalingFunctions -\> {None, "Log"}  x-axis linear, y-axis log ScalingFunctions -\> None           identity (default) ScalingFunctions -\> Automatic      identity (default) When a log scale is active, tick labels show original data-space values (e.g. 1, 2, 5, 10, 20, …) at decade-based positions. Non-positive values on a log-scaled axis are suppressed (mapped to -1e30 in world space).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= ScalingFunctions -> {"Log", None}
Out[1]= ScalingFunctions -> {"Log", None}
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[ScalingFunctions], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[Plot[Exp[x], {x, 1, 10}, ScalingFunctions -> "Log"]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `ScalingFunctions` is an inert option keyword for `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot` and `StreamPlot`: it applies a coordinate transform to one or both axes — "Log", "Log10", "Log2" or "Reverse" for both axes, or a {xform, yform} pair to select each axis independently; None or Automatic is the identity. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `ScalingFunctions -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `ScalingFunctions -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot` and `StreamPlot`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `ScalingFunctions` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ScalingFunctions` is an inert option keyword for `Plot`, `ListPlot`, `DensityPlot`, `ContourPlot`, `VectorPlot` and `StreamPlot` — a coordinate transform applied to one or both axes ("Log", "Log10", "Log2", "Reverse", or a per-axis pair). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ScalingFunctions -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
