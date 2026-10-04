---
source: src/graphics/graphics_init.c
---
**Definition.** `BarSpacing` is an inert option keyword for `BarChart` and `Histogram`: it sets the gap between bars as a fraction of bar width (default 0.2; 0 gives touching bars, 1 all gap). It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `BarSpacing -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `BarSpacing -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `BarChart` and `Histogram`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `BarSpacing` itself, which carries no logic of its own.
