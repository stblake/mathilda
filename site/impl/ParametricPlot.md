---
source: src/graphics/parametricplot.c
---
**Algorithm.** `builtin_parametricplot` is `HoldAll`. `is_iterator` tests each
trailing `{var, min, max}` spec; one iterator is the **curve** form, two is the
**filled region** form. The body is any expression that must evaluate to a
2-element list of finite reals — a literal `{fx, fy}` or a computed form such as
`r {Cos[t], Sin[t]}`. Two evaluation paths live in `ParamEvalCtx`: a compiled
fast path (`param_ctx_compile`, all-or-nothing, firing only for a *literal*
2-element `List` whose two coordinates each compile with `autocompile_new`) and
an interpreter fallback (`param_eval_at` → `eval_body_xy`, which binds the
iterator via an OwnValue and checks the result is a finite real pair). The
one-iterator form (`build_param_curve` → `param_sample`) runs its **own**
adaptive bisection (`param_subdivide`), the same three-probe refinement as `Plot`
but measuring **Euclidean chord deviation normalised by the bounding-box
diagonal** (neither axis is privileged for a parametric curve), deliberately kept
at the same tolerance `PARAM_FLAT_TOL = 0.0006`; invalid or region-rejected
samples break the polyline into separate `Line[...]` runs (one coloured two-point
`Line` per segment under `ColorFunction`). The two-iterator form
(`build_param_region`) samples a **uniform `n x n`** `(t, r)` grid, maps each pair
to `(x, y)`, and emits one filled `Polygon[{p00, p10, p11, p01}]` per cell whose
four corners are all valid — no adaptivity — with `Mesh -> All` overlaying the
grid edges as `Line`s. A `List`-of-`List`s body is the multi-curve form, each
sub-body drawn in a `palette_color` with non-colour directives given their own
`List` scope. The result is an inert `Graphics[prims, opts...]`.

**Data structures.** `ParamEvalCtx` (iterator vars, body, `acx`/`acy` compiled
coordinates, `RegionFunction`); `ParamPt`/`ParamBuf` growing sample buffers; the
option bundle from `split_options_param`.

**Complexity / limits.** `PlotPoints` initial samples (default **50** for the
curve, **75** for the region) plus refinement to `MaxRecursion = 6` for the
curve; the region grid is uniform with no refinement and drops any cell with an
invalid corner (no partial-triangle handling). Default `AspectRatio -> 1`,
`Axes -> True`, `PlotStyle` colour `RGBColor[0.2, 0.4, 0.8]`.
