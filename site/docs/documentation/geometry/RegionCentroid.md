# RegionCentroid

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RegionCentroid[Polygon[{{x1, y1}, ...}]] gives the centroid {cx, cy} of a simple 2D polygon with nonzero area. Exact coordinates give an exact result; Real coordinates give a machine-precision result.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[1]= {1/3, 1/3}

In[2]:= RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[2]= {2, 7/5}
```

### Applications (2)

The centroid of a unit right triangle

```mathematica
In[3]:= RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[3]= {1/3, 1/3}
```

Area-weighted, and exact

```mathematica
In[4]:= RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[4]= {2, 7/5}
```

## Implementation notes

**Algorithm.** `builtin_region_centroid` unwraps `Polygon[pts]`, reads and dedups the vertices,
and declines (returns `NULL`, leaving the call unevaluated) for fewer than 3 distinct
vertices. The polygon centroid uses the standard shoelace-weighted formula
`c_x = Sum_i (x_i + x_{i+1}) cross_i / (6A)`, `c_y = Sum_i (y_i + y_{i+1}) cross_i / (6A)`,
where `cross_i = x_i y_{i+1} - x_{i+1} y_i` and `2A` is the signed shoelace area:

- **exact path** — `geom_signed_area2_q` gives `2A`; if it is zero the polygon is degenerate
  and `RegionCentroid` **declines** (a documented deviation from WL, which returns the
  lower-dimensional measure centroid). Otherwise the two weighted sums are accumulated in GMP
  `mpq_t`, divided by `6A` (built as `3 * 2A`), and each component is canonicalised with
  `make_rational_mpz`, so `RegionCentroid[Polygon[{{0,0},{1,0},{0,1}}]]` is `{1/3, 1/3}`.
- **machine path** — the same formula in `double`, declining when the `double` signed area is
  exactly `0.0`.

**Data structures.** A `GeomPoints` with parallel `mpq_t`/`double` coordinate arrays; a block
of seven `mpq_t` temporaries (`a2`, `sx`, `sy`, `cr`, `t1`, `t2`, `six`) initialised once per
call on the exact path. The result is a two-element `List` built by `geom_list2`.

**Complexity / limits.** `O(n)` in the vertex count, plus GMP cost on the exact path. Simple
polygons only; `Protected`, not `Listable`. Declines on zero area, and (machine path) on a
coordinate too large for a `double`. Gates use exact comparisons — the zero-area test compares
to exactly `0.0`, no epsilon.

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)

## Notes & additional examples

### Notes

`RegionCentroid[Polygon[{{x1, y1}, …}]]` is the area centroid `{cx, cy}` of a simple 2D
polygon, from the shoelace-weighted formula. Exact coordinates give an exact `{cx, cy}`; any
`Real` coordinate gives a machine result (contagion).

It requires a polygon of **nonzero area**. A degenerate (zero-area) polygon makes
`RegionCentroid` **decline** and stay unevaluated — a deliberate deviation from the Wolfram
Language, which returns the lower-dimensional measure centroid (e.g. a segment midpoint). It
also declines on fewer than three distinct vertices. The zero-area gate is an exact comparison
(no epsilon). `Protected`, not `Listable`.
