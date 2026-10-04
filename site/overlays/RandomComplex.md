### Worked examples

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[]  (* uniform in the unit square [0,1]+[0,1] I *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[2 + 2 I]  (* rectangle with corners 0 and 2 + 2 I *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[{0, 1 + I}, 3]  (* three points in a rectangle *)
```

```mathematica
In[1]:= SeedRandom[1]; RandomComplex[{-1 - I, 1 + I}, {2, 2}]  (* a 2x2 array of points *)
```

### Notes

A complex range is an **axis-aligned rectangle** in the plane, not a disk: the real
and imaginary parts are drawn independently and uniformly. A bare `z` gives the
rectangle with corners `0` and `z`; `{z1, z2}` gives the rectangle spanned by `z1`
and `z2`.

Each draw is two uniforms from the same stream `RandomReal` uses, so `SeedRandom`
makes the sequence reproducible — the examples above are stable across runs.
`WorkingPrecision -> d` draws the parts at `d` digits through MPFR. A list or array
form (`RandomComplex[range, n]`, `RandomComplex[range, {n1, ...}]`) returns the
corresponding nested list of `Complex` values.
