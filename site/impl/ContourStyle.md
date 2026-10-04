---
source: src/graphics/graphics_init.c
---
**Definition.** `ContourStyle` is an inert option keyword for `ContourPlot`: it gives the style directive(s) for the contour lines; a single directive applies to all levels, a list cycles through them, Automatic colours by height and None suppresses the lines. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `ContourStyle -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `ContourStyle -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `ContourPlot`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `ContourStyle` itself, which carries no logic of its own.
