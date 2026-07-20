# BarChart

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
BarChart[{v1, v2, ..., vn}, opts...]
    Draws a vertical bar chart: n bars at x = 1..n with heights v1..vn.
BarChart[{{v1,...}, {w1,...}, ...}, opts...]
    Multiple grouped datasets, each in a distinct palette colour.

    Options:
      ChartStyle    color/style list cycling through bars (default: palette)
      ChartLabels   list of x-axis tick labels
      BarSpacing    gap fraction of bar width (default 0.2)
      Standard Graphics options (Axes, AspectRatio, Frame, PlotRange,
      PlotLabel, Background, ImageSize, …) pass through.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
