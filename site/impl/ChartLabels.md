---
source: src/graphics/graphics_init.c
---
**Definition.** `ChartLabels` is an **option name** for `BarChart` giving the labels drawn
below the bars on the x-axis. It is a bare option symbol, not a function: `graphics_init`
(`src/graphics/graphics_init.c`) stamps it `Protected` and sets its docstring, but it has
no builtin and no DownValues. The `BarChart` renderer reads it from the call's option list.

**Representation.** `ChartLabels` stays an inert, `Protected` `EXPR_SYMBOL`, appearing only
on the left of a rule; `FullForm[ChartLabels -> {"Q1", "Q2"}]` is
`Rule[ChartLabels, List["Q1", "Q2"]]`. Its value is a list of label expressions, one per
bar, in the bar order.

**Usage & limits.** Supplied as `BarChart[{3, 5, 2}, ChartLabels -> {"a", "b", "c"}]`, it
is consumed when the chart's x-axis is drawn and otherwise does nothing; it is a `BarChart`
option, not a general `Graphics` directive. It affects only the rendered labelling, never
the plotted data.
