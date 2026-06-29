# $PlotResample

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
$PlotResample[var, {f...}, {plotPoints, maxRecursion, maxPlotPoints, mesh, regionFunction, exclusions, colorFunction, colorFunctionScaling, filling, fillingStyle}]
    Internal Plot metadata used by the renderer to re-sample curves at the current zoom. Not intended for direct use.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
