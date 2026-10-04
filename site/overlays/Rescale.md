### Worked examples

```mathematica
In[1]:= Rescale[5, {0, 10}]  (* the midpoint maps to 1/2 *)
```

```mathematica
In[1]:= Rescale[Range[5]]  (* one argument uses the data's own Min and Max *)
```

```mathematica
In[1]:= Rescale[{1, 2, 3, 4}, {0, 10}, {0, 100}]  (* rescale into an explicit target range *)
```

### Notes

`Rescale[x, {min, max}]` maps `x` linearly so that `min -> 0` and `max -> 1`:
`(x - min)/(max - min)`. The three-argument form
`Rescale[x, {min, max}, {y0, y1}]` maps into an arbitrary target interval, and
the one-argument `Rescale[list]` uses `{Min[list], Max[list]}` as the source
range.

`Rescale` threads over lists at every level and keeps input exact — `Rescale[Range[5]]`
is `{0, 1/4, 1/2, 3/4, 1}`, not floats — because it rewrites its call into a
single `Plus`/`Times` expression over the whole argument and lets the `Listable`
arithmetic heads do the work, which also routes a packed array through the fast
vectorised kernels.
