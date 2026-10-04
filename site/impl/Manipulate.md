---
source: src/graphics/manipulate.c
---
**Algorithm.** `builtin_manipulate` is `HoldAll` and `Protected`. Like `Animate`
it returns C `NULL` (leaving `Manipulate[...]` symbolic, body unevaluated) when
`MATHILDA_NO_WINDOW` is set, so doc-generation never opens a window. Otherwise
each control spec (args 1+, each evaluated first) is parsed by
`parse_control_spec` into one of up to `MAX_CTRLS = 8` controls:
`parse_var_part` reads a bare symbol `u` or a `{u, u0}` default pair; a single
remaining nonempty `List` makes a **discrete** button set (`MCTRL_DISCRETE`, a
packed/visible `NDArray` value list is first materialised with
`ndarray_to_nested_list`, capped at `MAX_DISCRETE_VALUES = 32`), while two or
three further numeric args make a **continuous** range (`MCTRL_RANGE`,
numericised through `N[]`, with an optional step `du` that enables snapping). It
then calls `graphics_manipulate` and returns `Null`. Per frame the binding is
cleaner than `Animate`'s: for each control it `iter_spec_shadow`s any prior value,
installs an OwnValue (`expr_new_real(value)` for a range, a copy of the selected
value for a discrete set), evaluates `expr_copy(body)`, then restores in reverse
order. One control row is drawn per variable (drag sliders for ranges, clickable
button sets for discrete controls); a footer **Reset** button restores every
control's default and the 3-D camera. 2-D/3-D dispatch is identical to `Animate`
(`graphics_render_in_region` / `graphics3d_render_in_region` with an orbit
`cam3d`). There is no playback transport. The Raylib loop blocks the REPL until
`Esc` closes it.

**Data structures.** `ManipCtrl` (a `MCTRL_RANGE`/`MCTRL_DISCRETE` kind, the var
`Expr*`, range fields `vmin`/`vmax`/`step`/`value`/`default_value`, discrete
fields `values`/`n_values`/`selected_idx`/`default_idx`); loop state
`drag_ctrl`, `Expr* frame_expr`, `Graphics3DEmbedState* cam3d`. Window is
800×500, `MAX_CTRLS = 8`, `MAX_DISCRETE_VALUES = 32`.

**Complexity / limits.** One full body evaluation **every frame** (`SetTargetFPS(60)`),
unconditionally — the "whenever a control changes" wording in the docstring is
looser than the code. Headless/no-Raylib behaviour matches `Animate` (a one-line
message, still returns `Null`). At most 8 controls and 32 discrete values.
