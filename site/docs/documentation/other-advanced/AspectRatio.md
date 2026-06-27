# AspectRatio

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
AspectRatio
    is an option for Graphics and Plot that specifies the ratio of
    height to width of the rendered plot.

AspectRatio -> Automatic sets the ratio from the actual coordinate
values (true geometry); AspectRatio -> Full stretches the graphics to
fill the enclosing region; AspectRatio -> a uses the explicit
height-to-width ratio a. Plot defaults to 1/GoldenRatio.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
