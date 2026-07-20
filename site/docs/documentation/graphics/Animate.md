# Animate

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Animate[expr, {t, tmin, tmax}, opts...]
    Opens an interactive animation window, evaluating expr at each
    frame with t bound to the current parameter value. Returns Null
    once the window is closed. expr is typically a Graphics[...] or
    Plot[...] call that depends on t.

    Options:
      AnimationDirection    Forward (default) | Backward |
                            ForwardBackward | BackwardForward
      AnimationRate         parameter units per second (real > 0)
      AnimationRepetitions  integer or Infinity (default Infinity)
      AnimationRunning      True (default) | False (start paused)
      AppearanceElements    All (default) | None |
                            {"PlayPauseButton", "ProgressSlider",
                             "StepLeftButton", "StepRightButton",
                             "DirectionButton",
                             "FasterSlowerButtons", "ResetButton"}
      DefaultDuration       seconds for one full pass (default 1.0)
      ControlPlacement      Bottom (default) | Top
      RefreshRate           target display FPS (default 60)

    Keyboard controls: Space (play/pause), Arrow keys (step),
      R (reset), Esc (close). Direction/speed buttons in the
      control bar are clickable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
