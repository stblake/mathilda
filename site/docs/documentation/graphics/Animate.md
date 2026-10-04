# Animate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Animate[expr, {t, tmin, tmax}, opts...]`**

Opens an interactive animation window, evaluating expr at each frame with t bound to the current parameter value. Returns Null once the window is closed. expr is typically a Graphics\[...\] or Plot\[...\] call that depends on t. Options: AnimationDirection    Forward (default) | Backward | ForwardBackward | BackwardForward AnimationRate         parameter units per second (real \> 0) AnimationRepetitions  integer or Infinity (default Infinity) AnimationRunning      True (default) | False (start paused) AppearanceElements    All (default) | None | {"PlayPauseButton", "ProgressSlider", "StepLeftButton", "StepRightButton", "DirectionButton", "FasterSlowerButtons", "ResetButton"} DefaultDuration       seconds for one full pass (default 1.0) ControlPlacement      Bottom (default) | Top RefreshRate           target display FPS (default 60) Keyboard controls: Space (play/pause), Arrow keys (step), R (reset), Esc (close). Direction/speed buttons in the control bar are clickable.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Animate[Plot[Sin[x + t], {x, 0, 2 Pi}], {t, 0, 2 Pi}]
Out[1]= Animate[Plot[Sin[x + t], {x, 0, 2 Pi}], {t, 0, 2 Pi}]

In[2]:= Animate[Plot3D[Sin[x + t] Cos[y], {x, -3, 3}, {y, -3, 3}], {t, 0, 2 Pi}]
Out[2]= Animate[Plot3D[Sin[x + t] Cos[y], {x, -3, 3}, {y, -3, 3}], {t, 0, 2 Pi}]
```

### Options (4)

```mathematica
In[3]:= Animate[Graphics[Disk[{t, 0}, 0.5], PlotRange -> {{0, 5}, {-1, 1}}], {t, 0, 5}, AnimationDirection -> ForwardBackward]
Out[3]= Animate[-Graphics-, {t, 0, 5}, AnimationDirection -> ForwardBackward]

In[4]:= Animate[ ParametricPlot[{Cos[u], Sin[u]}, {u, 0, t}], {t, 0.01, 2 Pi}, DefaultDuration -> 4, AnimationRunning -> False]
Out[4]= Animate[ParametricPlot[{Cos[u], Sin[u]}, {u, 0, t}], {t, 0.01, 2 Pi}, DefaultDuration -> 4, AnimationRunning -> False]

In[5]:= Animate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}, AnimationRepetitions -> 3, AnimationRate -> 2]
Out[5]= Animate[Plot[Sin[n x], {x, 0, 2 Pi}], {n, 1, 5}, AnimationRepetitions -> 3, AnimationRate -> 2]

In[6]:= Animate[Plot[Exp[-t x^2], {x, -3, 3}], {t, 0.1, 5}, ControlPlacement -> Top, AppearanceElements -> {"PlayPauseButton", "ProgressSlider"}]
Out[6]= Animate[Plot[Exp[-t x^2], {x, -3, 3}], {t, 0.1, 5}, ControlPlacement -> Top, AppearanceElements -> {"PlayPauseButton", "ProgressSlider"}]
```

### Applications (4)

```mathematica
In[7]:= Attributes[Animate]
Out[7]= {HoldAll, Protected}

In[8]:= Animate[Plot[Sin[x + t], {x, 0, 2 Pi}], {t, 0, 2 Pi}]
Out[8]= Animate[Plot[Sin[x + t], {x, 0, 2 Pi}], {t, 0, 2 Pi}]

In[9]:= Animate[Graphics[Disk[{t, 0}, 0.5], PlotRange -> {{0, 5}, {-1, 1}}], {t, 0, 5}, AnimationDirection -> ForwardBackward]
Out[9]= Animate[-Graphics-, {t, 0, 5}, AnimationDirection -> ForwardBackward]

In[10]:= Animate[Plot3D[Sin[x + t] Cos[y], {x, -3, 3}, {y, -3, 3}], {t, 0, 2 Pi}]
Out[10]= Animate[Plot3D[Sin[x + t] Cos[y], {x, -3, 3}, {y, -3, 3}], {t, 0, 2 Pi}]
```

## Algorithm

animate.c — Animate[expr, {t, tmin, tmax}, opts...]

Opens a Raylib window (when USE_GRAPHICS is compiled in) with a Manipulate-style control bar:

```text
  ┌─────────────────────────────────────────────────────────────┐
  │ a  │ 0 ─────────────●────────── 5  │  2.31               │
  ├─────────────────────────────────────────────────────────────┤
  │ [R][|<][▶][>|][dir]                              1x [-][+] │
  │ Space:play/pause  ←/→:step  R:reset  Esc:close            │
  └─────────────────────────────────────────────────────────────┘
```

Each positional iterator argument {var, min, max} adds a labeled

```text
slider row.  All iterators share a common animation phase so they
advance together.  Additional iterator arguments beyond the first
```

work as simultaneously animated variables (not static sliders).

Options supported:

```text
  AnimationDirection    Forward | Backward | ForwardBackward |
                        BackwardForward (default Forward)
  AnimationRate         parameter units per second of the FIRST iterator
                        (overrides DefaultDuration)
  AnimationRepetitions  integer or Infinity (default Infinity)
  AnimationRunning      True | False (default True)
  AppearanceElements    All | None | {"PlayPauseButton", ...}
  DefaultDuration       seconds for one full pass (default 1.0)
  ControlPlacement      Bottom | Top (default Bottom)
  RefreshRate           target display FPS (default 60)
```

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Show](../../graphics/Show/)

- Source: [`src/graphics/animate.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/animate.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`Animate` is `HoldAll` and `Protected`: the body `expr` is held and re-evaluated
every frame with the iterator variable bound to the current parameter value, a
single phase `phi` driving all iterators. It opens a Raylib window with playback
controls — `AnimationDirection`, `AnimationRate`, `AnimationRepetitions`,
`AnimationRunning`, `DefaultDuration`, `ControlPlacement`, `RefreshRate` — and the
event loop **blocks the REPL** until the window closes, whereupon it returns `Null`.
A `Graphics3D`/`Plot3D` body gets its own persistent orbit camera.

Under headless conditions (the `MATHILDA_NO_WINDOW` environment variable, as used by
documentation generation) `Animate` returns unevaluated — the body is never run and
no window opens — which is why the examples above display as the held
`Animate[...]` form here rather than opening a window.
