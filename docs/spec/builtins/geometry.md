# Computational geometry

Core 2D computational-geometry heads over `Polygon` and point sets, following
Wolfram Language 12 semantics (every example below is locked to a live Wolfram
kernel reference output, 2026-08-27; ticket GEO-1). Implemented in
`src/geometry.c`; registered by `geometry_init()`.

All five heads have attributes `{Protected, ReadProtected}` (WL-faithful) and
are deliberately **not** `Listable` — they are structural, not element-wise.

Two computation paths, chosen per call:

- **Exact** — every coordinate is `Integer`/`BigInt`/`Rational`: all arithmetic
  in GMP rationals; results are exact (`Integer`/`Rational`, or symbolic
  radical sums for `Perimeter`).
- **Machine** — any `Real`/`MPFR` coordinate switches the whole computation to
  machine doubles (WL's contagion semantics). A visible `NDArray` points
  argument (`Polygon[NDArray[...]]`, `ConvexHullRegion[NDArray[...]]`, or an
  `NDArray` query point) is read directly and is always machine.

Anything out of scope declines and **stays unevaluated**: symbolic
coordinates, non-2D points, `Polygon` with holes (the two-argument component
form), other region heads (`Disk`, `Rectangle`, ...), wrong arity.

## Area

`Area[Polygon[{{x1, y1}, ...}]]` — shoelace area of a simple 2D polygon.

```
Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]         (* 1/4 *)
Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]   (* 10, concave *)
Area[Polygon[{{0, 0}, {1.5, 0}, {1.5, 1}, {0, 1}}]] (* 1.5 *)
Area[Polygon[{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}}]]   (* 1: explicit closing vertex OK *)
Area[Polygon[{{0, 0}, {1, 0}}]]                     (* Undefined: < 3 distinct vertices *)
```

## Perimeter

`Perimeter[Polygon[{{x1, y1}, ...}]]` — sum of edge lengths including the
closing edge. Exact input builds `Sqrt` terms the evaluator canonicalizes.

```
Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]     (* 2 + Sqrt[2] *)
Perimeter[Polygon[{{0, 0}, {1/2, 0}, {0, 1}}]]   (* 3/2 + 1/2 Sqrt[5] *)
Perimeter[Polygon[{{0, 0}, {3., 0}, {3., 4.}}]]  (* 12.0 *)
Perimeter[Polygon[{{0, 0}, {1, 0}}]]             (* Undefined *)
```

## RegionCentroid

`RegionCentroid[Polygon[{{x1, y1}, ...}]]` — area centroid `{cx, cy}` of a
simple polygon with nonzero area.

```
RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]              (* {1/3, 1/3} *)
RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]  (* {2, 7/5} *)
```

## RegionMember

`RegionMember[Polygon[{{x1, y1}, ...}], {x, y}]` — `True` if the point is
inside or **on the boundary** of the simple polygon (vertices included),
`False` otherwise. Exact input uses exact sign tests; boundary decisions on
exact input are exact.

```
RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]     (* True  *)
RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {2, 1}]     (* True: on edge *)
RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {3, 1}]     (* False *)
RegionMember[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}], {2, 3}]  (* False: in the notch *)
```

## ConvexHullRegion

`ConvexHullRegion[{{x1, y1}, ...}]` — convex hull of a 2D point set (Andrew's
monotone chain; duplicate, interior, and collinear-middle points removed).
Result by hull dimension: `Polygon` (vertices counterclockwise, starting from
the lexicographic minimum — WL's own vertex order), `Line` for collinear
input, `Point` for a single point.

```
ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]
  (* Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}] *)
ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]   (* Line[{{0, 0}, {3, 3}}] *)
ConvexHullRegion[{{1, 2}}]                           (* Point[{1, 2}] *)
Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]  (* 4 *)
```

## When these heads decline

Declining leaves the expression unevaluated, which is always preferable to a
confident wrong number. In addition to symbolic coordinates and non-2D input:

- **A coordinate too large for a double, asked about on a machine path.**
  `RegionMember[Polygon[{{0,0},{10^400,0},{0,1}}], {1., 1.}]` declines: the
  exact coordinate mirrors to infinity, every cross product becomes NaN, and a
  NaN reads as "collinear" — which produced `True`, the wrong answer, before
  this check existed. The fully exact spelling (`{1, 1}`) still answers `False`.
  Likewise `ConvexHullRegion[{{0,0},{10^400,0},{0,1},{1.,1.}}]` declines rather
  than returning a degenerate `Line` containing a non-re-parseable `inf.0`.
- **Fewer than 3 distinct vertices**, for `RegionCentroid` — and `Undefined`
  for `Area`/`Perimeter`. Consecutive repeats (including the wrap-around pair)
  are collapsed first, so `Area[Polygon[{{0,0},{0,0},{1,1}}]]` is `Undefined`,
  not `0`, and `Perimeter[Polygon[{{0,0},{1,0},{1,0}}]]` is `Undefined`, not a
  there-and-back `2`.
- **A degenerate (zero-area) polygon**, for `RegionMember` and `RegionCentroid`.

## Exactness across representations

| Input form | Result | Why |
|---|---|---|
| nested `List` of exact numbers | exact | exact GMP path |
| **packed** integer list (e.g. from `Table`) | **exact** | packed buffers materialise before these heads run, which is what preserves exactness — the reason they are `EXEMPT` rather than `AWARE` in `tools/check_packed_aware.py` |
| visible `NDArray[...]` | machine | `NDArray[{{0,0},...}]` stores and prints `0.0` — float64 by construction, so there is no exactness to keep. Matches `Total[NDArray[{1,2,3}]]` being `6.0` where `Total[{1,2,3}]` is `6` |
| any `Real` coordinate | machine | WL contagion |

## Documented deviations from Wolfram Language

- **Simple polygons only.** Self-intersecting vertex lists follow
  shoelace/even-odd semantics; WL computes the enclosed region. Not detected
  at runtime.
- **`RegionCentroid` declines on zero-area polygons** (WL returns the
  lower-dimensional measure centroid, e.g. the segment midpoint).
- **`ConvexHullRegion` returns bare `Polygon[{verts}]`** — WL 12+ attaches a
  cell-spec second argument (`Polygon[verts, {1, 2, ..., n}]`).
- **Machine-path comparisons are exact IEEE** — boundary membership requires
  the cross product to be exactly zero; the zero-area centroid gate compares
  to exactly `0.0`. No epsilon tolerances.
- Machine reals print in mathilda's form (`12.0`, `4.0`) where WL prints
  `12.`, `4.` — values identical.

## Tests

`tests/test_geometry.c` (ctest target `geometry_tests`) covers every example
above plus decline paths and a valgrind memory loop;
`tests/scripts/geometry_e2e.m` runs the same acceptance rows through the full
REPL pipeline (success marker: `geometry_e2e: ALL PASS` — callers must grep,
the process exit code is not a verdict).
