---
source: src/graphics/graphics_init.c
---
**Definition.** `ChartStyle` is an inert option keyword for `BarChart` and `Histogram`: it gives a colour, or a list of colours cycled through the bars; it defaults to the standard palette. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `ChartStyle -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `ChartStyle -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `BarChart` and `Histogram`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `ChartStyle` itself, which carries no logic of its own.
