# Perimeter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Perimeter[Polygon[{{x1, y1}, ...}]] gives the perimeter of a simple 2D polygon: the sum of its edge lengths, including the closing edge. Exact coordinates give an exact (possibly symbolic, e.g. 2 + Sqrt[2]) result; Real coordinates give a machine-precision result.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[1]= 2 + Sqrt[2]

In[2]:= Perimeter[Polygon[{{0, 0}, {1/2, 0}, {0, 1}}]]
Out[2]= 3/2 + 1/2 Sqrt[5]

In[3]:= Perimeter[Polygon[{{0, 0}, {3., 0}, {3., 4.}}]]
Out[3]= 12.0
```

### Applications (3)

An exact perimeter keeps the radical

```mathematica
In[4]:= Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[4]= 2 + Sqrt[2]
```

Rational coordinates, a canonicalised Sqrt sum

```mathematica
In[5]:= Perimeter[Polygon[{{0, 0}, {1/2, 0}, {0, 1}}]]
Out[5]= 3/2 + 1/2 Sqrt[5]
```

A real coordinate gives a machine length

```mathematica
In[6]:= Perimeter[Polygon[{{0, 0}, {3., 0}, {3., 4.}}]]
Out[6]= 12.0
```

## Implementation notes

**Algorithm.** `builtin_perimeter` unwraps `Polygon[pts]` with `geom_polygon_points`, reads
and dedups the vertices (`geom_read_points` + `geom_dedup_consecutive`, consecutive and
wrap-around repeats collapsed), and returns `Undefined` for fewer than 3 distinct vertices.
The perimeter is the sum of edge lengths including the closing edge,
`Sum_i Sqrt[(x_{i+1}-x_i)^2 + (y_{i+1}-y_i)^2]`:

- **exact path** — each squared edge length is computed in GMP (`mpq_sub`/`mpq_mul`/
  `mpq_add`), wrapped as a `Sqrt[...]` `Expr`, and the `Plus` of those terms is handed to the
  evaluator, so radical canonicalisation does the rest: `Perimeter[Polygon[{{0,0},{1,0},
  {0,1}}]]` becomes `2 + Sqrt[2]` and `Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]]` becomes
  `3/2 + 1/2 Sqrt[5]`.
- **machine path** — a running `double` sum of `hypot(dx, dy)`, returned as a `Real`.

The exact/machine choice is made once by `geom_read_points` from the coordinate types (WL
contagion: any `Real` coordinate makes the whole result machine).

**Data structures.** A `GeomPoints` with parallel `mpq_t qx/qy` and `double dx/dy`. The exact
path allocates one `Sqrt[...]` `Expr` per edge into a `terms` array, builds a single `Plus`,
and lets `eval_and_free` simplify it (an OOM mid-build frees the partial terms and declines).

**Complexity / limits.** `O(n)` edges, plus the evaluator's radical-simplification cost on the
exact path and GMP arithmetic on each squared length. Simple polygons only; `Protected`, not
`Listable`. A coordinate that overflows a `double` on a machine path makes `Perimeter`
decline.

**Attributes:** `Protected`.

## References

**See also:** [Sqrt](../../arithmetic/Sqrt/)

- Source: [`src/geometry.c`](https://github.com/stblake/mathilda/blob/main/src/geometry.c)
- Specification: [`docs/spec/builtins/geometry.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/geometry.md)
- Tests: [`tests/test_geometry.c`](https://github.com/stblake/mathilda/blob/main/tests/test_geometry.c)

## Notes & additional examples

### Notes

`Perimeter[Polygon[{{x1, y1}, …}]]` is the sum of the edge lengths, including the closing edge
from the last vertex back to the first.

On exact (`Integer`/`Rational`) coordinates each edge length is a `Sqrt` of an exact squared
distance, and the sum is handed to the evaluator, so radical canonicalisation produces the
Wolfram-shaped exact answer (`2 + Sqrt[2]`, `3/2 + 1/2 Sqrt[5]`, …). Any `Real` coordinate
makes the whole result a machine `Real`, summed with `hypot`.

Like `Area`, consecutive and wrap-around duplicate vertices are collapsed first, and a polygon
with fewer than three distinct vertices is `Undefined`. Simple polygons only; `Protected`, not
`Listable`.
