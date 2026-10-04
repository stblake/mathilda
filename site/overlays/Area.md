### Worked examples

```mathematica
In[1]:= Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]  (* exact rational coordinates give an exact area *)
```

```mathematica
In[1]:= Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]  (* a concave polygon, by the shoelace sum *)
```

```mathematica
In[1]:= Area[Polygon[{{0, 0}, {1.5, 0}, {1.5, 1}, {0, 1}}]]  (* any real coordinate makes the result machine *)
```

```mathematica
In[1]:= Area[Polygon[{{0, 0}, {1, 0}}]]  (* fewer than three distinct vertices is Undefined *)
```

### Notes

`Area[Polygon[{{x1, y1}, …}]]` is the shoelace area of a simple 2D polygon.

The coordinate types choose the path, once per call. All-`Integer`/`Rational` coordinates
compute in exact GMP rationals, so the answer is an exact `Integer` or `Rational`. Any `Real`
coordinate switches the whole computation to machine doubles (the Wolfram Language's
contagion rule), giving a `Real`.

Consecutive duplicate vertices — including an explicit closing vertex that repeats the first —
are collapsed before counting, so a closed vertex list works and a polygon left with fewer
than three distinct vertices is `Undefined` rather than a confident `0`.

Simple (non-self-intersecting) polygons only: a self-intersecting vertex list follows
shoelace/even-odd semantics, which differ from the enclosed-region model, and that is not
detected at runtime. `Area` is `Protected` and not `Listable`. A visible `NDArray` of
`int64` points keeps the exact path; a float `NDArray` is machine.
