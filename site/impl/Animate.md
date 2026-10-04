---
source: src/graphics/animate.c
---
**Algorithm.** `builtin_animate` is `HoldAll` and `Protected`, so the body `expr`
and the iterator spec `{t, tmin, tmax}` arrive unheld only after this builtin
chooses to evaluate them. It first checks `getenv("MATHILDA_NO_WINDOW")`: if set,
it returns C `NULL` immediately — the head stays the symbolic `Animate[...]`, the
body is never evaluated, and no window opens (this is what lets doc-generation
re-run without popping windows). Otherwise it collects up to `MAX_ITERS = 8`
iterator specs (each evaluated and numericised through `N[]`, bounds swapped if
inverted), parses the playback options with `parse_animate_opts`, calls
`graphics_animate`, and returns the `Null` symbol once the window closes. Per
frame, `graphics_animate` does **binding, not substitution**: for each iterator it
installs a global OwnValue `var = <real>` via `symtab_add_own_value`, calls
`evaluate(expr_copy(body))`, then clears the symbol. A single phase `phi ∈ [0,1]`
drives all iterators (`t = tmin + phi*span`); `AnimationDirection` governs loop
vs ping-pong at the endpoints, honouring `AnimationRepetitions` (`Infinity` = −1).
The evaluated frame is dispatched by head: `Graphics` →
`graphics_render_in_region`, `Graphics3D` → `graphics3d_render_in_region` with a
persistent orbit camera, anything else → text via `expr_to_string`. The window
runs a synchronous Raylib event loop (`while (!WindowShouldClose())`) on the
calling thread, so it **blocks the REPL** until closed.

**Data structures.** `AnimIter {const char* var_sym; double tmin, tmax, t}`;
`AnimateOpts` (direction enum, rate, repetitions, running, duration, control
placement, refresh rate, and five AppearanceElements flags); per-run loop state
`phi`, `going_fwd`, `running`, `speed_mult` (clamped `[0.125, 8]`), the current
`Expr* frame_expr`, and a `Graphics3DEmbedState* cam3d`. Window is 800×500.

**Complexity / limits.** One full body evaluation **every rendered frame** (at
`RefreshRate`), unconditionally — not only when the parameter changes. Headless
behaviour: `MATHILDA_NO_WINDOW` leaves the call unevaluated, while
`MATHILDA_NO_GRAPHICS_WINDOW` and a `USE_GRAPHICS=0` build print a one-line
message and return `Null`. `"ProgressSlider"`/`"AnimationSlider"` are accepted
`AppearanceElements` names but currently have no effect; the per-iterator slider
rows are always drawn.
