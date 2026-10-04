# LightGreen

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightGreen`**

The named colour RGBColor\[0.88, 1, 0.88\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A pale pastel -- its own constant, not a tint of Green

```mathematica
In[1]:= LightGreen
Out[1]= RGBColor[0.88, 1, 0.88]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[LightGreen]
Out[2]= RGBColor[0.88, 1, 0.88]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{LightGreen, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `LightGreen` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightGreen", rgb_color(0.88, 1, 0.88))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[0.88, 1, 0.88]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightGreen` evaluates to `RGBColor[0.88, 1, 0.88]`, so
`FullForm[LightGreen]` is `RGBColor[0.88, 1, 0.88]` and `Head[LightGreen]` is `RGBColor`.
The exact `1` green component prints as an integer while the `0.88` reds print as reals
(`color_component`). It is a separate constant from `Green`, one of the pale `Light*`
pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightGreen, Disk[]}]`, `PlotStyle -> LightGreen`, a light `ArrayPlot` cell. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightGreen` is the named colour `RGBColor[0.88, 1, 0.88]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightGreen]` is `RGBColor[0.88, 1, 0.88]` and `Head[LightGreen]` is `RGBColor`.

The exact `1` green component prints as an integer while the `0.88` reds print as reals.
It is a separate constant from `Green`, one of the pale `Light*` pastel family. It carries
no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
