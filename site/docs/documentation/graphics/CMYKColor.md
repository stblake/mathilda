# CMYKColor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CMYKColor[c, m, y, k]`**

represents a color in the CMYK (cyan, magenta, yellow, black) space.

**`CMYKColor[c, m, y, k, a] specifies opacity a; CMYKColor[c, m, y] takes`**

<details>
<summary>Notes</summary>

k = 0. The list forms CMYKColor\[{c, m, y, k}\] and CMYKColor\[{c, m, y, k, a}\] are also accepted. Components and opacity outside \[0,1\] are clipped. A style directive: sets the colour of subsequent graphics primitives, converted to RGB as r=(1-c)(1-k), g=(1-m)(1-k), b=(1-y)(1-k).

</details>

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

A cyan disk

```mathematica
In[1]:= Graphics[{CMYKColor[1, 0, 0, 0], Disk[]}]
Out[1]= -Graphics-
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
