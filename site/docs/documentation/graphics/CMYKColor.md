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

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

A cyan disk

```mathematica
In[1]:= Graphics[{CMYKColor[1, 0, 0, 0], Disk[]}]
Out[1]= -Graphics-
```

### Applications (4)

```mathematica
In[2]:= CMYKColor[1, 0, 0, 0]
Out[2]= CMYKColor[1, 0, 0, 0]

In[3]:= Graphics[{CMYKColor[1, 0, 0, 0], Disk[]}]
Out[3]= -Graphics-

In[4]:= MemberQ[Attributes[CMYKColor], Protected]
Out[4]= True

In[5]:= Head[CMYKColor[0.2, 0.4, 0.1, 0.3]]
Out[5]= CMYKColor
```

## Implementation notes

**Algorithm.** `CMYKColor` is an inert, `Protected` **style directive** for the
subtractive (print) colour model — it has no builtin and no OwnValue (only the
attribute bit is set in `graphics_init.c`; its docstring lives centrally in
`info.c`), so it stays symbolic and is interpreted only at render time, exactly
like `RGBColor`/`GrayLevel`/`Hue`. The forms are `CMYKColor[c, m, y, k]`,
`CMYKColor[c, m, y]` (`k = 0`), `CMYKColor[c, m, y, k, a]` (with opacity), and the
list forms `CMYKColor[{c, m, y, k}]` / `CMYKColor[{c, m, y, k, a}]`. At render
time `rgba_from_cmyk`/`cmyk_to_rgb` (`render.c`) clip every component to `[0,1]`
(`clip01`) and convert with `w = 1 - k`, `r = (1 - c) w`, `g = (1 - m) w`,
`b = (1 - y) w`. It is recognised as a colour directive by the interned-pointer
test `h == SYM_CMYKColor` in `resolve_color` wherever `RGBColor`/`GrayLevel`/`Hue`
are — as a primitive-list directive, a `PlotStyle`/`Background`/`FrameStyle`
value, or a `ColorFunction` return value. The **same** conversion is duplicated in
`graphics_json.c` for the JSON/notebook-frontend path, so the CLI Raylib renderer
and the notebook agree.

**Data structures.** No dedicated CMYK struct: it is converted straight to an
`RGBA8 {unsigned char r, g, b, a}` at parse/render time (and thence to Raylib's
`Color` via `to_raylib`), or serialised to an `rgba(...)` string in the JSON
path.

**Complexity / limits.** `O(1)` per colour. Components and opacity outside
`[0,1]` are clipped rather than rejected. Because it is inert, `CMYKColor[...]`
never evaluates to an `RGBColor` at the language level — the conversion happens
only inside the two renderers.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/render.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/render.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

### Notes

`CMYKColor` is an inert, `Protected` style directive for the subtractive (print)
colour model. It has no builtin and no own-values, so it never evaluates to an
`RGBColor` at the language level — it stays symbolic and is converted only at render
time, by `r = (1 - c)(1 - k)`, `g = (1 - m)(1 - k)`, `b = (1 - y)(1 - k)`. The forms
are `CMYKColor[c, m, y, k]`, `CMYKColor[c, m, y]` (`k = 0`), `CMYKColor[c, m, y, k,
a]` (with opacity), and the list forms `CMYKColor[{c, m, y, k}]` /
`CMYKColor[{c, m, y, k, a}]`; components and opacity outside `[0, 1]` are clipped.

It is recognised as a colour wherever `RGBColor`/`GrayLevel`/`Hue` are — as a
primitive-list directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a
`ColorFunction` return value — and the same conversion is implemented in both the
CLI Raylib renderer and the notebook JSON frontend.
