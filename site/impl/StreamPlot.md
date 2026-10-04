---
references:
  - "B. Jobard and W. Lefebvre, *Creating Evenly-Spaced Streamlines of Arbitrary Density*, Proc. 8th Eurographics Workshop on Visualization in Scientific Computing (1997) 43-55."
source: src/graphics/streamplot.c
---
**Algorithm.** `builtin_streamplot` is `HoldAll`. It traces **evenly-spaced**
streamlines of the field `{vx, vy}` by the Jobard–Lefebvre scheme. Each line is
grown in **both directions** from a seed (`grow_streamline`/`integrate_dir`) by
`rk4_step_unit`, a 4th-order Runge–Kutta step of the **normalised (unit) field**
`û = v/|v|` with signed step `s = ±h` — integrating the unit field gives uniform
arc-length point spacing independent of local speed. With `ext = min(dx, dy)` and
`np = StreamPoints` (default **25**) the data-space parameters are: target line
separation `d_sep = ext/np`; RK4 step `h = d_sep/8`; proximity-stop radius²
`d_test2 = 0.25*d_sep^2` (a line halts ~`0.5*d_sep` from an existing line);
closed-orbit return radius² `loop2 = (1.6 h)^2` after a minimum arc `6 h`. Even
spacing is enforced by a uniform spatial hash `SGrid` with **cell size exactly
`d_sep`**: every accepted streamline point is registered, and
`sgrid_within2(x, y, r2)` searches the query cell's 3×3 neighbourhood (correct
because `r2 <= cell^2`), giving `O(1)` proximity tests. Seeds are drawn from a
candidate grid twice as fine as `d_sep` (`2*np` per axis), sorted by distance from
the domain centre so placement grows outward. A line stops at the domain
boundary, a `RegionFunction` rejection, a critical point (`|v| <= crit`),
proximity to an existing line, a closed-orbit return, or the step cap
(`max_arc = StreamScale*diag` or `4*diag`, `max_steps` capped at 20000). Default
rendering (`emit_dashed_stream`) walks arc length alternating drawn `dash_len`
and blank `gap_len`, each dash a curve-following `Line` capped by a filled
`Polygon` arrowhead triangle (`chevron_tip`) — a bare triangle, not `Arrow`, to
avoid a straight chord shaft. `StreamAnimate -> True` instead emits one
`AnimatedStreamline[...]` per line plus periodic `Arrow` chevrons, whose particle
dots flow in the interactive renderer. Colour priority is
`StreamColorFunction`/`ColorFunction` at the line midpoint, else a Viridis
speed ramp, else a single `StreamStyle` directive. The result is an inert
`Graphics[prims, opts...]`.

**Data structures.** `SGrid` (the `d_sep`-cell spatial hash: `nx*ny` buckets of
`Point2` with per-cell `cnt`/`cap`); `SeedCand {x, y, d2}` sorted candidate array;
`PBuf` growing per-half point buffers; `all_streams`/`all_lengths`/`all_speed`
parallel arrays. Integration runs in data space; finished polylines map to world
space only for rendering.

**Complexity / limits.** `(2*np)^2` seed candidates (default 2500), `O(1)`
proximity per point via the hash, per-line steps capped at 20000. Returns `NULL`
if no streamline is produced. Defaults `Frame -> True`, `Axes -> False`,
`AspectRatio -> 1`.
