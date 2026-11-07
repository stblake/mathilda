# Computational geometry

Core 2D computational-geometry heads over `Polygon` and point sets, following
the reference semantics (every example below is locked to a live reference
kernel reference output, 2026-08-27; ticket GEO-1). Implemented in
`src/geometry.c`; registered by `geometry_init()`.

All five heads have attributes `{Protected}` (reference-faithful) and
are deliberately **not** `Listable` — they are structural, not element-wise.

Two computation paths, chosen per call:

- **Exact** — every coordinate is `Integer`/`BigInt`/`Rational`: all arithmetic
  in GMP rationals; results are exact (`Integer`/`Rational`, or symbolic
  radical sums for `Perimeter`).
- **Machine** — any `Real`/`MPFR` coordinate switches the whole computation to
  machine doubles (the reference's contagion semantics). A visible `NDArray` points
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

```
In[1]:= Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]
Out[1]= 1/4

In[2]:= Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[2]= 10

In[3]:= Area[Polygon[{{0, 0}, {1.5, 0}, {1.5, 1}, {0, 1}}]]
Out[3]= 1.5

In[4]:= Area[Polygon[{{0, 0}, {1, 0}}]]
Out[4]= Undefined
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

```
In[1]:= Perimeter[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[1]= 2 + Sqrt[2]

In[2]:= Perimeter[Polygon[{{0, 0}, {1/2, 0}, {0, 1}}]]
Out[2]= 3/2 + 1/2 Sqrt[5]

In[3]:= Perimeter[Polygon[{{0, 0}, {3., 0}, {3., 4.}}]]
Out[3]= 12.0
```

## RegionCentroid

`RegionCentroid[Polygon[{{x1, y1}, ...}]]` — area centroid `{cx, cy}` of a
simple polygon with nonzero area.

```
RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]              (* {1/3, 1/3} *)
RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]  (* {2, 7/5} *)
```

```
In[1]:= RegionCentroid[Polygon[{{0, 0}, {1, 0}, {0, 1}}]]
Out[1]= {1/3, 1/3}

In[2]:= RegionCentroid[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
Out[2]= {2, 7/5}
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

```
In[1]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {1, 1}]
Out[1]= True

In[2]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {2, 1}]
Out[2]= True

In[3]:= RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {3, 1}]
Out[3]= False
```

## ConvexHullRegion

`ConvexHullRegion[{{x1, y1}, ...}]` — convex hull of a 2D point set (Andrew's
monotone chain; duplicate, interior, and collinear-middle points removed).
Result by hull dimension: `Polygon` (vertices counterclockwise, starting from
the lexicographic minimum — the reference's own vertex order), `Line` for collinear
input, `Point` for a single point.

```
ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]
  (* Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}] *)
ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]   (* Line[{{0, 0}, {3, 3}}] *)
ConvexHullRegion[{{1, 2}}]                           (* Point[{1, 2}] *)
Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]  (* 4 *)
```

```
In[1]:= ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]
Out[1]= Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]

In[2]:= ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]
Out[2]= Line[{{0, 0}, {3, 3}}]

In[3]:= ConvexHullRegion[{{1, 2}}]
Out[3]= Point[{1, 2}]

In[4]:= Area[ConvexHullRegion[{{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}}]]
Out[4]= 4
```

## ArcLength

`ArcLength[reg]` / `ArcLength[curve, {t, a, b}]` (src/geometry_curves.c).

- `Line[{p1, ..., pn}]`: the sum of segment lengths, any dimension
  (`ArcLength[Line[{{0,0},{3,4},{3,5}}]] = 6`); a list of point lists adds the
  pieces. Exact coordinates give the exact sum of square roots; machine
  coordinates sum scaled `hypot` over the doubles directly (a 10^5-point contour
  takes ~50 ms and builds no expression); a packed point list is read through
  its nested-List form so an int64 buffer stays exact. `WorkingPrecision -> p`
  gives `N[sum, p]`.
- `Circle[c, r]`: `2 Pi r`; an ellipse `Circle[c, {a, b}]` gives
  `4 b EllipticE[1 - a^2/b^2]` with `b` the larger semi-axis
  (`8 EllipticE[3/4]` for `{1, 2}`); an arc `Circle[c, r, {t1, t2}]` gives
  `r (t2 - t1)`.
- A parametrised curve `{x1[t], ..., xn[t]}` integrates `|c'(t)|`; a single
  function `f[t]` its graph, `Sqrt[1 + f'(t)^2]`. The squared speed is
  simplified first (`{Cos[t], Sin[t]}` over `{t, 0, 2 Pi}` is `2 Pi`), then
  integrated exactly (`ArcLength[{t, t^2}, {t, 0, 1}] = (2 Sqrt[5] + ArcSinh[2])/4`).
  An inexact bound or `WorkingPrecision -> p` integrates by `NIntegrate`.
- **Deliberate difference:** when `Integrate` finds no closed form over a numeric
  range, `ArcLength` returns the `NIntegrate` value (machine precision); the
  reference returns the inert `Integrate[...]`, but `N` of an inert `Integrate`
  is not a number in Mathilda, so that form would be a dead end. Elliptic cases
  the reference closes (`ArcLength[Sin[t], {t, 0, Pi}] = 2 Sqrt[2] EllipticE[1/2]`)
  therefore come back as `3.82019778902771`. A symbolic range keeps the inert
  integral. (The exact attempt costs whatever `Integrate` takes to give up, a few
  seconds for `Sqrt[1 + 3 Cos[t]^2]`.)
- `Point[...]`, `Disk[...]`, `Annulus[...]`: `Undefined`.
- Attributes: `Protected`.

## RegionMeasure

`RegionMeasure[reg]`: `Area` of a `Polygon`, `Disk` or `Annulus`; `ArcLength`
of a `Line` or `Circle`. Attributes: `Protected`.

## Annulus and Disk

`Annulus[{x, y}, {rin, rout}]` is the closed annulus, `Annulus[c, {rin, rout},
{t1, t2}]` the annular sector; `Annulus[]` is `Annulus[{0, 0}, {1/2, 1}]` (and
`Annulus[c]` uses the same radii). Attributes: `Protected`.

- `Area` / `RegionMeasure`: `Pi (rout^2 - rin^2)` (`3 Pi/4` for `Annulus[]`;
  `Pi (-r1^2 + r2^2)` symbolically); a sector `(t2 - t1)/2 (rout^2 - rin^2)`.
  Numeric radii given out of order measure as the ordered annulus
  (`Area[Annulus[{0,0},{2,1}]] = 3 Pi`), as in the reference.
- `Perimeter`: `2 Pi (rin + rout)`; a sector `(t2 - t1)(rin + rout) + 2 (rout - rin)`.
- `RegionCentroid`: the centre (full annulus).
- `RegionMember[ann, p]`: `rin <= |p - c| <= rout` (closed), plus the angle of
  `p - c` within `[t1, t2]` for a sector; a symbolic point gives the reference's
  condition `Element[x | y, Reals] && 1/4 <= x^2 + y^2 <= 1`.
- `Disk[c, r]` / `Disk[c, {a, b}]` / `Disk[]`: `Area` / `RegionMeasure`
  `Pi r^2` / `Pi a b` / `Pi`; `Perimeter` the boundary circle's `ArcLength`;
  `RegionCentroid` the centre.
- Not implemented: `RegionBounds`, `RegionDimension`, the one-argument
  `RegionMember[reg]` (`RegionMemberFunction`), and use as an `NDSolve` /
  `NDEigenvalues` domain.

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
| any `Real` coordinate | machine | reference contagion |

## Documented deviations from the reference

- **Simple polygons only.** Self-intersecting vertex lists follow
  shoelace/even-odd semantics; the reference computes the enclosed region. Not detected
  at runtime.
- **`RegionCentroid` declines on zero-area polygons** (the reference returns the
  lower-dimensional measure centroid, e.g. the segment midpoint).
- **`ConvexHullRegion` returns bare `Polygon[{verts}]`** — the reference attaches a
  cell-spec second argument (`Polygon[verts, {1, 2, ..., n}]`).
- **Machine-path comparisons are exact IEEE** — boundary membership requires
  the cross product to be exactly zero; the zero-area centroid gate compares
  to exactly `0.0`. No epsilon tolerances.
- Machine reals print in mathilda's form (`12.0`, `4.0`) where the reference prints
  `12.`, `4.` — values identical.

## Tests

`tests/test_geometry.c` (ctest target `geometry_tests`) covers every example
above plus decline paths and a valgrind memory loop;
`tests/scripts/geometry_e2e.m` runs the same acceptance rows through the full
REPL pipeline (success marker: `geometry_e2e: ALL PASS` — callers must grep,
the process exit code is not a verdict).
