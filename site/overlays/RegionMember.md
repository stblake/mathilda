### Worked examples

```mathematica
In[1]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]  (* an interior point *)
```

```mathematica
In[1]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {2, 1}]  (* a boundary point counts as a member *)
```

```mathematica
In[1]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {3, 1}]  (* an exterior point *)
```

```mathematica
In[1]:= RegionMember[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}], {2, 3}]  (* a point in a concave notch is outside *)
```

### Notes

`RegionMember[Polygon[{{x1, y1}, …}], {x, y}]` gives `True` when the point is inside **or on
the boundary** of a simple 2D polygon, and `False` otherwise. The test is a crossing-number
ray cast with an explicit on-edge check, so vertices and edge points are members.

Exact coordinates and an exact query point use exact GMP sign tests, so boundary decisions are
exact. If either side is a `Real`, the machine path runs with exact IEEE comparisons (a cross
product of exactly `0.0` is on-edge — no epsilon tolerance). An exact polygon whose
coordinates do not fit a `double`, asked about with a machine point, declines rather than
answer from NaN cross products.

It declines on a degenerate (zero-area) polygon or fewer than three distinct vertices.
`Protected`, not `Listable`.
