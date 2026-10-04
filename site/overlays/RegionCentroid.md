### Worked examples

```mathematica
In[1]:= RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]  (* the centroid of a unit right triangle *)
```

```mathematica
In[1]:= RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]  (* area-weighted, and exact *)
```

### Notes

`RegionCentroid[Polygon[{{x1, y1}, …}]]` is the area centroid `{cx, cy}` of a simple 2D
polygon, from the shoelace-weighted formula. Exact coordinates give an exact `{cx, cy}`; any
`Real` coordinate gives a machine result (contagion).

It requires a polygon of **nonzero area**. A degenerate (zero-area) polygon makes
`RegionCentroid` **decline** and stay unevaluated — a deliberate deviation from the Wolfram
Language, which returns the lower-dimensional measure centroid (e.g. a segment midpoint). It
also declines on fewer than three distinct vertices. The zero-area gate is an exact comparison
(no epsilon). `Protected`, not `Listable`.
