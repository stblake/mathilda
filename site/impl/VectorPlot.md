---
source: src/graphics/vectorplot.c
---
**Algorithm.** `builtin_vectorplot` is `HoldAll`; `args[0]` must be a 2-element
list `{vx, vy}`. The field is compiled all-or-nothing by `vp_compile` into two
`AutoCompiled` kernels `f(x, y)` (`VpField`), with `vp_eval` falling back to the
interpreter per sample. It samples an **`N x N` grid** (`N = VectorPoints`,
default **15**) in scaled world space, inverting each point to data space
(`scale_invert`) before evaluation and masking it with `RegionFunction`; when a
`ScalingFunctions` transform is active the arrow direction is corrected by a
finite-difference Jacobian. Arrow sizing is **screen-normalised** so arrows stay
legible across mixed-scale axes: directions are divided by the per-axis world
units, and `half_len_screen = 0.5/(N-1)` (half a grid cell) sets the base length.
`VectorScale` selects the length rule — `Automatic` gives a fixed length,
`None` makes length proportional to magnitude (`len ∝ speed/speed_max`), a real
`f` scales relative to the grid spacing. Each arrow is **centred** on its grid
point (`tail = grid - d`, `head = grid + d`) and emitted as
`Arrow[{{tail}, {head}}]` preceded by its colour directive; `vp_color` defaults to
the `default_ramp_rgb` (Viridis) ramp keyed to the normalised speed (a user
function is tried at arity 3 `f[vx, vy, speed]` then arity 1 `f[speed]`). The
result is an inert `Graphics[prims, opts..., PlotRange -> {...}]`.

**Data structures.** `VpField` (the two compiled kernels); a local
`Arrow {wx, wy, vx, vy}` array of `N*N` world-space entries; the option bundle
from `split_vector_options`. The primitive list is a leading `VectorStyle` (or the
default `Thickness[0.005]`) then a colour directive and `Arrow` per grid point.

**Complexity / limits.** `O(N^2)` field evaluations and arrows (default 225). The
grid is always square (no per-axis point count), and only `Arrow` primitives are
emitted. Defaults `Frame -> True`, `Axes -> False`, `AspectRatio -> 1`,
`ColorFunctionScaling -> True`.
