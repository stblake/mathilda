### Worked examples

```mathematica
In[1]:= ListInterpolation[{1, 4, 9, 16}]  (* values on the integer grid 1, 2, 3, 4 *)
```

```mathematica
In[1]:= ListInterpolation[{1, 4, 9, 16}][2.5]  (* evaluate the interpolant between grid points *)
```

```mathematica
In[1]:= ListInterpolation[{1, 4, 9, 16}, {{0, 3}}]  (* place the grid equally spaced on [0, 3] *)
```

```mathematica
In[1]:= ListInterpolation[{1, 4, 9, 16}, InterpolationOrder -> 1][2.5]  (* piecewise-linear: halfway between 4 and 9 *)
```

```mathematica
In[1]:= ListInterpolation[{{1, 2}, {3, 4}}][1.5, 1.5]  (* a 2-D array interpolates in both directions *)
```

### Notes

`ListInterpolation[array]` builds an `InterpolatingFunction` from an array of
*values* taken to lie on a regular grid at integer positions `1, 2, ...` in each
direction — the value-only companion to `Interpolation`, which takes
`{abscissa, value}` pairs. The nesting depth of `array` is the number of
dimensions. Call the returned object like a function to sample it.

A second argument places the grid: `{{xmin, xmax}, ...}` spaces it equally over an
interval per dimension, or explicit position lists give the grid lines directly.
`InterpolationOrder -> n` sets the degree (default 3, so a longer run of points is
a cubic fit; `1` is piecewise-linear), and `Method` and `PeriodicInterpolation`
select the scheme.
