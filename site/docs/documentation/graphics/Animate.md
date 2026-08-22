# Animate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Animate[expr, {t, tmin, tmax}, opts...]`**

Opens an interactive animation window, evaluating expr at each frame with t bound to the current parameter value. Returns Null once the window is closed. expr is typically a Graphics\[...\] or Plot\[...\] call that depends on t. Options: AnimationDirection    Forward (default) | Backward | ForwardBackward | BackwardForward AnimationRate         parameter units per second (real \> 0) AnimationRepetitions  integer or Infinity (default Infinity) AnimationRunning      True (default) | False (start paused) AppearanceElements    All (default) | None | {"PlayPauseButton", "ProgressSlider", "StepLeftButton", "StepRightButton", "DirectionButton", "FasterSlowerButtons", "ResetButton"} DefaultDuration       seconds for one full pass (default 1.0) ControlPlacement      Bottom (default) | Top RefreshRate           target display FPS (default 60) Keyboard controls: Space (play/pause), Arrow keys (step), R (reset), Esc (close). Direction/speed buttons in the control bar are clickable.

## Examples (6)

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Show](../../graphics/Show/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
