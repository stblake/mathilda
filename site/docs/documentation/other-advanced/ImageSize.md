# ImageSize

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ImageSize
    is an option for Graphics and Plot that specifies the overall size
    of the image to display.

ImageSize -> w sets the width to w pixels, with the height following
from AspectRatio; ImageSize -> {w, h} fixes both the width and height
in pixels. The default is a width of 800.
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
