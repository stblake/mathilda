---
source: src/graphics/graphics_init.c
---
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
