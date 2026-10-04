### Worked examples

```mathematica
In[1]:= ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]  (* interior and collinear points are dropped *)
```

```mathematica
In[1]:= ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]  (* collinear input collapses to a Line *)
```

```mathematica
In[1]:= ConvexHullRegion[{{1, 2}}]  (* a single point gives a Point *)
```

```mathematica
In[1]:= Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]  (* composes with the other region heads *)
```

### Notes

`ConvexHullRegion[{{x1, y1}, …}]` is the convex hull of a 2D point set, computed by Andrew's
monotone chain. Duplicate, interior, and collinear-middle points are removed, and the vertices
come out counterclockwise starting from the lexicographic minimum — the Wolfram Language's own
vertex order.

The result head follows the hull's dimension: a `Polygon` for a genuine 2D hull, a `Line` for
collinear input, and a `Point` for a single distinct point. It returns a bare
`Polygon[{verts}]` — unlike Wolfram 12+, which attaches a cell-spec second argument. Integer/
rational input is processed in exact arithmetic (so the sign tests are exact); a `Real`
coordinate switches to the machine path. A point coordinate too large for a `double` makes the
machine path decline rather than emit a degenerate region. `Protected`, not `Listable`.
