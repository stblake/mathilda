# Area

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Area[Polygon[{{x1, y1}, ...}]] gives the area of a simple 2D polygon. Exact (Integer/Rational) coordinates give an exact result; any Real coordinate gives a machine-precision result. A polygon with fewer than 3 distinct vertices has Undefined area.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]
Out[1]= 1/4

In[2]:= Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[2]= 10

In[3]:= Area[Polygon[{{0, 0}, {1.5, 0}, {1.5, 1}, {0, 1}}]]
Out[3]= 1.5

In[4]:= Area[Polygon[{{0, 0}, {1, 0}}]]
Out[4]= Undefined
```

### Applications (4)

Exact rational coordinates give an exact area

```mathematica
In[5]:= Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]
Out[5]= 1/4
```

A concave polygon, by the shoelace sum

```mathematica
In[6]:= Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[6]= 10
```

Any real coordinate makes the result machine

```mathematica
In[7]:= Area[Polygon[{{0, 0}, {1.5, 0}, {1.5, 1}, {0, 1}}]]
Out[7]= 1.5
```

Fewer than three distinct vertices is Undefined

```mathematica
In[8]:= Area[Polygon[{{0, 0}, {1, 0}}]]
Out[8]= Undefined
```

## Implementation notes

**Algorithm.** `builtin_area` unwraps `Polygon[pts]` (exactly one part) via
`geom_polygon_points`, reads the vertices with `geom_read_points`, and collapses consecutive
duplicate vertices — including the wrap-around pair — with `geom_dedup_consecutive`, so an
explicitly closed vertex list and interior repeats are handled. A polygon with fewer than 3
distinct vertices is degenerate and returns `Undefined`. The area is the shoelace sum
`|Sum_i (x_i y_{i+1} - x_{i+1} y_i)| / 2`:

- **exact path** (every coordinate `Integer`/`BigInt`/`Rational`) — `geom_signed_area2_q`
  accumulates twice the signed area in a GMP `mpq_t`, which is then made positive, halved
  (`mpq_div_2exp`), and canonicalised to an `Integer` or `Rational` by `make_rational_mpz`. So
  `Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]]` is exactly `1/4`.
- **machine path** (any `Real`/`MPFR` coordinate, WL's contagion) — `geom_signed_area2_d`
  sums in `double` and the result is `fabs(area2)/2` as a `Real`.

`geom_read_points` classifies the whole vertex list in one pass (exact until a `Real` forces
machine; decline on a symbolic or complex coordinate), so the two paths never mix.

**Data structures.** A `GeomPoints` holding parallel `mpq_t qx/qy` (exact) and `double
dx/dy` (always mirrored), with an `exact` flag and a `d_finite` flag. A visible `NDArray` of
points is read directly: an `int64`-typed array keeps the exact path (so it agrees with the
same points as a nested `List`), any other dtype takes the machine path.

**Complexity / limits.** `O(n)` in the number of vertices, plus GMP cost on the exact path.
Simple (non-self-intersecting) polygons only — a self-intersecting vertex list follows
shoelace/even-odd semantics, not WL's enclosed-region model, and is not detected at runtime.
`Area` is `Protected` and deliberately **not** `Listable` (it is structural, not
element-wise); a coordinate too large for a `double` on a machine path makes it decline rather
than compute on an infinity.

**Attributes:** `Protected`.

## References

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)

## Notes & additional examples

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
