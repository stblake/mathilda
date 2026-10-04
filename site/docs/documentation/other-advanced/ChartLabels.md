# ChartLabels

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ChartLabels`**

BarChart option: list of label expressions drawn below each bar on the x-axis.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A bare BarChart option symbol

```mathematica
In[1]:= Head[ChartLabels]
Out[1]= Symbol
```

Inert and Protected -- no builtin behaviour

```mathematica
In[2]:= Attributes[ChartLabels]
Out[2]= {Protected}
```

It only sits on the left of a rule

```mathematica
In[3]:= FullForm[ChartLabels -> {"Q1", "Q2"}]
Out[3]= Rule[ChartLabels, List["Q1", "Q2"]]
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ChartLabels` is an option for `BarChart` giving the labels drawn below the bars on the
x-axis. It is an inert, `Protected` option-name symbol — no builtin, no DownValues — so it
does nothing on its own; the `BarChart` renderer reads it from the call's option list.

Its value is a list of label expressions, one per bar, in bar order, as in
`BarChart[{3, 5, 2}, ChartLabels -> {"a", "b", "c"}]`. It affects only the rendered
labelling, never the data, and is specific to `BarChart` rather than a general `Graphics`
directive.
