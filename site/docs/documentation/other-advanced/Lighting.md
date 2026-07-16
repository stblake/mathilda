# Lighting

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Lighting
    A Graphics3D / Plot3D / ParametricPlot3D option controlling surface shading. Lighting -> Automatic (default): per-face Lambertian (flat) shading with a fixed directional light; ambient 0.3, diffuse 0.7. Lighting -> None (or False): disables shading and draws surfaces in their raw PlotStyle/ColorFunction color.
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
