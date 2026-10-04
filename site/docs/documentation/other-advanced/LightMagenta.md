# LightMagenta

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightMagenta`**

The named colour RGBColor\[1, 0.9, 1\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

The named constant evaluates to its RGBColor literal

```mathematica
In[1]:= LightMagenta
Out[1]= RGBColor[1, 0.9, 1]
```

Its full form is a plain RGBColor triple

```mathematica
In[2]:= FullForm[LightMagenta]
Out[2]= RGBColor[1, 0.9, 1]
```

A protected style constant

```mathematica
In[3]:= MemberQ[Attributes[LightMagenta], Protected]
Out[3]= True
```

Used as a directive, styling the primitive that follows it

```mathematica
In[4]:= Graphics[{LightMagenta, Disk[]}]
Out[4]= -Graphics-
```

## Implementation notes

**Definition.** `LightMagenta` is a named colour constant standing for `RGBColor[1, 0.9, 1]` — a pale magenta tint
in the additive (screen) RGB model. It is **not** an inert head: `register_color`
(`graphics_init.c`) installs it as an **OwnValue**, so `LightMagenta` *evaluates* to its
`RGBColor[1, 0.9, 1]` literal, and the `Protected` attribute is set on the same symbol. The
docstring is kept centrally in `info.c`, mirroring `CMYKColor` and the other named
colours.

**Representation.** The OwnValue is built by `rgb_color(...)`, which runs each
channel through `color_component`: a component exactly `0` or `1` is emitted as a
plain `Integer` (matching Mathematica's `FullForm`), any other value as a `Real`.
Because `LightMagenta` reduces to a real colour literal before it ever reaches a renderer,
it resolves anywhere `RGBColor`/`GrayLevel`/`Hue` do — a primitive-list directive,
a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction` return.

**Usage & limits.** A style directive only: it sets the colour of the graphics
primitives that follow it and carries no numeric meaning elsewhere — `D`, `N` and
arithmetic see the `RGBColor[1, 0.9, 1]` triple, not a scalar. The named colours are defined
together as one table in `graphics_init.c`; `LightMagenta` is a single row of it.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightMagenta` is a named colour constant for `RGBColor[1, 0.9, 1]` (a pale magenta tint). Unlike `CMYKColor`, which is
an inert head, `LightMagenta` is an **OwnValue**: it *evaluates* to its `RGBColor[...]`
literal, which is why it resolves anywhere a real colour would — as a primitive-list
directive, a `PlotStyle`/`Background`/`FrameStyle` value, or a `ColorFunction`
result. A channel of exactly `0` or `1` prints as an integer; other channels print
as reals. It is `Protected`, and is one of the named colours defined together in
`graphics_init.c`.
