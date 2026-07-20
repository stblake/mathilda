# ScalingFunctions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ScalingFunctions
    Option for Plot, ListPlot, DensityPlot, ContourPlot, VectorPlot,
    StreamPlot: applies a coordinate transform to one or both axes.

    Forms:
      ScalingFunctions -> "Log"         both axes: natural log
      ScalingFunctions -> "Log10"        both axes: log base 10
      ScalingFunctions -> "Log2"         both axes: log base 2
      ScalingFunctions -> "Reverse"      both axes: mirror (negate)
      ScalingFunctions -> {"Log", None}  x-axis log, y-axis linear
      ScalingFunctions -> {None, "Log"}  x-axis linear, y-axis log
      ScalingFunctions -> None           identity (default)
      ScalingFunctions -> Automatic      identity (default)

    When a log scale is active, tick labels show original data-space
    values (e.g. 1, 2, 5, 10, 20, …) at decade-based positions.
    Non-positive values on a log-scaled axis are suppressed (mapped to
    -1e30 in world space).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
