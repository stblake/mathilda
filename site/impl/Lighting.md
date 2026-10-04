---
source: src/graphics/graphics_init.c
---
**Definition.** `Lighting` is an inert option keyword for `Graphics3D`, `Plot3D` and `ParametricPlot3D`: it controls surface shading for 3D graphics; Automatic (the default) is per-face Lambertian shading with a fixed directional light, while None (or False) disables shading and draws the raw PlotStyle/ColorFunction colour. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `Lighting -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `Lighting -> value` rule by the owning plot
builtin while it assembles the `Graphics3D[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `Graphics3D`, `Plot3D` and `ParametricPlot3D`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `Lighting` itself, which carries no logic of its own.
