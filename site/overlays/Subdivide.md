### Worked examples

```mathematica
In[1]:= Subdivide[5]  (* 0 to 1 in five equal parts, exact rationals *)
```

```mathematica
In[1]:= Subdivide[0, 10, 5]  (* five parts spanning 0 to 10 *)
```

```mathematica
In[1]:= Subdivide[0, 1.0, 4]  (* a machine-real endpoint makes every point a real *)
```

```mathematica
In[1]:= Subdivide[5, 5, 2]  (* a degenerate interval repeats the endpoint *)
```

### Notes

`Subdivide[n]` gives `n + 1` equally spaced points from `0` to `1`;
`Subdivide[max, n]` spans `0` to `max`, and `Subdivide[min, max, n]` spans `min`
to `max`. `n` counts the *parts*, so there are always `n + 1` points. A
descending interval (`max < min`) just produces descending points, with no
special case.

Results are exact when the endpoints are: `Subdivide[5]` is
`{0, 1/5, 2/5, 3/5, 4/5, 1}`, and whole points print as integers alongside the
rationals. A machine-real endpoint makes the whole result machine reals (and, at
scale, a packed array). The endpoints are copied from the input, never recomputed,
and each interior point is derived directly from its index, so nothing drifts.
