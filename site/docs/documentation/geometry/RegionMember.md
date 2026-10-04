# RegionMember

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegionMember[Polygon[{{x1, y1}, ...}], {x, y}] gives True if the point lies inside or on the boundary of the simple 2D polygon, and False otherwise.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]
Out[1]= True

In[2]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {2, 1}]
Out[2]= True

In[3]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {3, 1}]
Out[3]= False
```

### Applications (4)

An interior point

```mathematica
In[4]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]
Out[4]= True
```

A boundary point counts as a member

```mathematica
In[5]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {2, 1}]
Out[5]= True
```

An exterior point

```mathematica
In[6]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {3, 1}]
Out[6]= False
```

A point in a concave notch is outside

```mathematica
In[7]:= RegionMember[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}], {2, 3}]
Out[7]= False
```

## Implementation notes

**Algorithm.** `builtin_region_member` takes `RegionMember[Polygon[pts], {x, y}]`, reads and
dedups the polygon vertices, declines for fewer than 3 distinct vertices, and reads the query
point with `geom_read_query`. It is a boundary-inclusive point-in-polygon test by the
crossing-number (ray-casting) method, with an explicit on-edge check:

- **exact path** (`geom_member_q`, used when both polygon and query are exact) — for each edge
  it takes the sign of the cross product `(v_{i+1}-v_i) × (p-v_i)` with `geom_cross_sign_q`
  (exact GMP). A zero cross with the point inside the edge's bounding box means the point lies
  **on** the boundary → immediately `True`. Otherwise an upward/downward crossing of the
  horizontal ray toggles the inside flag, with the crossing direction matched to the cross
  sign so a vertex is counted once. The vertices themselves count as members.
- **machine path** (`geom_member_d`) — the identical logic in `double`, using **exact IEEE
  comparisons** (a cross product of exactly `0.0` is on-edge; no epsilon tolerance).

It returns the symbol `True` or `False`. An exact polygon whose coordinates do not fit a
`double`, queried with a machine point, **declines** rather than compute NaN cross products
that would read as "collinear" and answer wrongly.

**Data structures.** A `GeomPoints` with parallel `mpq_t`/`double` coordinates and a
`d_finite` flag; a `GeomCrossTmp` block of six `mpq_t` scratch values initialised once and
reused across the whole edge walk (the exact cross-sign predicate sits in an `O(n)` loop, so
per-call `mpq_init`/`mpq_clear` would dominate it).

**Complexity / limits.** `O(n)` edges per query. Simple polygons only; `Protected`, not
`Listable`. Declines on a zero-area/degenerate polygon and on the oversized-coordinate machine
case above. Boundary membership is exact on exact input.

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)

## Notes & additional examples

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
