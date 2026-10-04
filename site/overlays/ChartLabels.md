### Worked examples

```mathematica
In[1]:= Head[ChartLabels]  (* a bare BarChart option symbol *)
```

```mathematica
In[1]:= Attributes[ChartLabels]  (* inert and Protected -- no builtin behaviour *)
```

```mathematica
In[1]:= FullForm[ChartLabels -> {"Q1", "Q2"}]  (* it only sits on the left of a rule *)
```

### Notes

`ChartLabels` is an option for `BarChart` giving the labels drawn below the bars on the
x-axis. It is an inert, `Protected` option-name symbol — no builtin, no DownValues — so it
does nothing on its own; the `BarChart` renderer reads it from the call's option list.

Its value is a list of label expressions, one per bar, in bar order, as in
`BarChart[{3, 5, 2}, ChartLabels -> {"a", "b", "c"}]`. It affects only the rendered
labelling, never the data, and is specific to `BarChart` rather than a general `Graphics`
directive.
