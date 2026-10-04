---
source: src/graphics/graphics_init.c
---
**Algorithm.** `PolarAxes` is not a plotter — it is an inert, `Protected` option
keyword for `PolarPlot`, registered in `graphics_init.c` with only an attribute
bit and a docstring (there is no builtin function and no OwnValue). A grep of the
whole tree finds `PolarAxes` referenced nowhere outside `graphics_init.c`: it is
not one of the names `split_options_param` consumes, so when passed as
`PolarAxes -> True` it falls into the generic option pass-through and rides
verbatim as a trailing `Rule` onto the resulting `Graphics[...]` object. The
renderer does not yet draw the requested polar grid (radial circles plus angular
degree/radian labels), so **Cartesian axes are drawn instead** — `PolarPlot`
injects the usual `Axes -> True`. The symbol therefore stays symbolic in every
context and has no evaluation rules of its own.

**Data structures.** None — it is a bare interned symbol carrying
`ATTR_PROTECTED` and a docstring; its value, when supplied, is stored as an
ordinary `Rule[PolarAxes, True]` option on the graphics object.

**Complexity / limits.** None. The feature is deliberately a documented
placeholder: the option is recognised and preserved so existing `PolarPlot`
calls that set it keep working, but the polar grid overlay is unimplemented and
the plot renders with Cartesian axes. This is recorded here rather than silently
dropped.
