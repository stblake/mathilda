### Worked examples

```mathematica
(* a single dataset: n bars at x = 1..n *)
In[1]:= BarChart[{3, 1, 4, 1, 5, 9, 2, 6}]
```

```mathematica
(* grouped datasets in distinct palette colours *)
In[1]:= BarChart[{{1, 3, 2}, {4, 2, 5}}, BarSpacing -> 0.3]
```

```mathematica
(* negative values draw downward from the axis *)
In[1]:= BarChart[{3, -1, 4, -1, 5}]
```

```mathematica
(* custom colours and category labels *)
In[1]:= BarChart[{2.5, 4.1, 3.3, 5.7}, ChartStyle -> {Red, Blue, Green, Orange}, ChartLabels -> {"Q1", "Q2", "Q3", "Q4"}]
```

```mathematica
(* one bar is three primitives plus two outline directives: n bars from n values *)
In[1]:= Length[Cases[BarChart[{3, 1, 4}], _Rectangle, Infinity]]
```

### Notes

`BarChart` is `Protected` but **not** `HoldAll`, so its data is evaluated and each
element coerced to a real (non-numeric entries are dropped). A flat list gives `n`
bars at `x = 1..n`; a list of lists gives grouped sub-bars, one `palette_color` per
dataset. Each bar spans `[min(v, 0), max(v, 0)]`, so negatives draw downward.

`BarSpacing` (default 0.2) is the gap as a fraction of the group width; `ChartStyle`
cycles explicit colours; `ChartLabels` attaches category labels below each bar. Each
bar is drawn as a filled `Rectangle` plus a thin outline. Defaults are
`Axes -> True`, `AspectRatio -> 0.618`.
