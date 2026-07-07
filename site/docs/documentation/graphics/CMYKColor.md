# CMYKColor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
CMYKColor[c, m, y, k]
    represents a color in the CMYK (cyan, magenta, yellow, black) space.
CMYKColor[c, m, y, k, a] specifies opacity a; CMYKColor[c, m, y] takes
k = 0. The list forms CMYKColor[{c, m, y, k}] and CMYKColor[{c, m, y, k,
a}] are also accepted. Components and opacity outside [0,1] are clipped.
A style directive: sets the colour of subsequent graphics primitives,
converted to RGB as r=(1-c)(1-k), g=(1-m)(1-k), b=(1-y)(1-k).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
