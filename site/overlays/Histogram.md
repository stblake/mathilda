### Worked examples

```mathematica
(* default binning by Sturges' rule *)
In[1]:= Histogram[{1, 2, 2, 3, 3, 3, 4, 4, 4, 4}]
```

```mathematica
(* an explicit bin count *)
In[1]:= Histogram[Range[100], 10]
```

```mathematica
(* a fixed bin width *)
In[1]:= Histogram[Range[100], {25}]
```

```mathematica
(* an explicit range and width: {min, max, step} *)
In[1]:= Histogram[Range[100], {0, 100, 25}]
```

```mathematica
(* 100 integers into 10 equal bins gives 10 populated bars *)
In[1]:= Length[Cases[Histogram[Range[100], 10], _Rectangle, Infinity]]
```

### Notes

`Histogram` lives alongside `BarChart` and shares its bar-drawing machinery; it is
`Protected` but **not** `HoldAll`. The bin specification is flexible:
`Histogram[data]` uses **Sturges' rule** (`ceil(log2 n) + 1` bins, clamped to
`[2, 50]`); `Histogram[data, k]` gives `k` bins; `Histogram[data, {step}]` fixes the
bin width; `Histogram[data, {min, max, step}]` fixes the range and width.

Each value is dropped into its bin and the bin counts become bar heights; empty
bins are skipped (no zero-height rectangle). Only equal-width frequency bins are
supported — there is no density/PDF normalisation. Defaults are `Axes -> True`,
`AspectRatio -> 0.618`.
