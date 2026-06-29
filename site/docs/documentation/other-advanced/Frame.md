# Frame

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Frame
    is an option for Graphics and Plot that specifies whether to draw a
    frame with ticks and labels around the plot.

Frame -> True boxes all four edges; Frame -> False (or None) draws no
frame; Frame -> {{left, right}, {bottom, top}} toggles each edge with
True or False. In Plot a frame takes the place of the default Axes.
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
