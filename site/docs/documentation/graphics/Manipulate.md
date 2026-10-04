# Manipulate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Manipulate[expr, {u, umin, umax}, ...]`**

Opens an interactive window with one control per variable and re-evaluates expr with each variable bound to its current value whenever a control changes. Returns Null once the window is closed. expr is typically a Graphics\[...\] or Plot\[...\] call that depends on the control variables. Unlike Animate, every control is independently user-driven from the first frame -- there is no animation phase or playback transport. Control specs (any number, one row each): {u, umin, umax}              continuous slider, default = umin {u, umin, umax, du}          continuous slider with step du {{u, u0}, umin, umax}        continuous slider, explicit default u0 {{u, u0}, umin, umax, du}    continuous slider, default + step {u, {v1, v2, ...}}           discrete button set, default = v1 {{u, u0}, {v1, v2, ...}}     discrete button set, explicit default u0 Click-drag a slider handle or click a discrete button to change its value; click Reset to restore every control to its default. Esc closes the window.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]
Out[1]= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]

In[2]:= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]
Out[2]= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]

In[3]:= Manipulate[Plot[a Sin[x] + b, {x, 0, 2 Pi}], {a, 0, 3}, {b, {-1, 0, 1}}]
Out[3]= Manipulate[Plot[a Sin[x] + b, {x, 0, 2 Pi}], {a, 0, 3}, {b, {-1, 0, 1}}]

In[4]:= Manipulate[Plot3D[Sin[x + n] Cos[y], {x, -3, 3}, {y, -3, 3}], {n, 0, 3}]
Out[4]= Manipulate[Plot3D[Sin[x + n] Cos[y], {x, -3, 3}, {y, -3, 3}], {n, 0, 3}]
```

### Options (1)

```mathematica
In[5]:= Manipulate[ Graphics[Disk[{0, 0}, r], PlotRange -> {{-5, 5}, {-5, 5}}], {{r, 2}, 0.5, 5, 0.25}]
Out[5]= Manipulate[-Graphics-, {{r, 2}, 0.5, 5, 0.25}]
```

### Applications (4)

```mathematica
In[6]:= Attributes[Manipulate]
Out[6]= {HoldAll, Protected}

In[7]:= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]
Out[7]= Manipulate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}]

In[8]:= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]
Out[8]= Manipulate[Plot[f, {x, -5, 5}], {f, {Sin[x], Cos[x], x^2}}]

In[9]:= Manipulate[Graphics[Disk[{0, 0}, r], PlotRange -> {{-5, 5}, {-5, 5}}], {{r, 2}, 0.5, 5, 0.25}]
Out[9]= Manipulate[-Graphics-, {{r, 2}, 0.5, 5, 0.25}]
```

## Algorithm

manipulate.c — Manipulate[expr, {u, umin, umax}, ...]

Opens a Raylib window (when USE_GRAPHICS is compiled in) with one control row per variable, stacked above a thin footer:

```text
  ┌─────────────────────────────────────────────────────────────┐
  │                                                               │
  │                        (rendered content)                    │
  │                                                               │
  ├─────────────────────────────────────────────────────────────┤
  │ u  │ 0 ─────────────●────────── 10              │  3.42     │
  │ f  │ [Sin[x]] [Cos[x]] [Tan[x]]                              │
  ├─────────────────────────────────────────────────────────────┤
  │ [Reset]  Esc: close                                          │
  └─────────────────────────────────────────────────────────────┘
```

Each control is independently user-driven from the first frame — there is no animation phase or playback transport (see Animate for that). A control is either:

```text
  - a continuous drag slider:  {u, umin, umax}, {u, umin, umax, du},
    {{u, u0}, umin, umax}, {{u, u0}, umin, umax, du}
  - a discrete button set:     {u, {v1, v2, ...}}, {{u, u0}, {v1, ...}}
```

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Animate](../../graphics/Animate/), [Show](../../graphics/Show/)

- Source: [`src/graphics/manipulate.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/manipulate.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_manipulate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_manipulate.c)

## Notes & additional examples

### Notes

`Manipulate` is `HoldAll` and `Protected`: one control row is drawn per variable and
the body is re-evaluated with each variable bound to its current value as a control
changes. Control specs cover continuous sliders (`{u, min, max}`, with an optional
step `du`, or a `{{u, u0}, ...}` explicit default) and discrete button sets
(`{u, {v1, v2, ...}}`). A footer **Reset** button restores every control's default.
Unlike `Animate` there is no playback transport — every control is user-driven.

The window blocks the REPL until closed, then returns `Null`; a `Graphics3D`/`Plot3D`
body gets its own orbit camera. Under the headless `MATHILDA_NO_WINDOW` environment
(used by documentation generation) `Manipulate` returns unevaluated, which is why the
examples above display as the held `Manipulate[...]` form rather than opening a
window.
