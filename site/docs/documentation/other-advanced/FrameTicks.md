# FrameTicks

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
FrameTicks
    is an option for Graphics and Plot that specifies the tick marks on
    the edges of a frame.

FrameTicks -> Automatic (the default) draws major and minor ticks with
labels on every drawn frame edge; FrameTicks -> None keeps the frame
box but draws no ticks; the {{left, right}, {bottom, top}} form selects
Automatic or None per edge.
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
