# LightYellow

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightYellow`**

The named colour RGBColor\[1, 1, 0.85\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A pale pastel -- distinct from Yellow

```mathematica
In[1]:= LightYellow
Out[1]= RGBColor[1, 1, 0.85]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[LightYellow]
Out[2]= RGBColor[1, 1, 0.85]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{LightYellow, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `LightYellow` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightYellow", rgb_color(1, 1, 0.85))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[1, 1, 0.85]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightYellow` evaluates to `RGBColor[1, 1, 0.85]`, so
`FullForm[LightYellow]` is `RGBColor[1, 1, 0.85]` and `Head[LightYellow]` is `RGBColor`.
The exact `1` red/green components print as integers and the `0.85` blue as a real
(`color_component`). It is a distinct constant from `Yellow` (`RGBColor[1, 1, 0]`), one of
the pale `Light*` pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightYellow, Disk[]}]`, `PlotStyle -> LightYellow`, a pale highlight fill. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightYellow` is the named colour `RGBColor[1, 1, 0.85]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightYellow]` is `RGBColor[1, 1, 0.85]` and `Head[LightYellow]` is `RGBColor`.

The exact `1` red/green components print as integers and the `0.85` blue as a real. It is
a distinct constant from `Yellow` (`RGBColor[1, 1, 0]`), one of the pale `Light*` pastel
family. It carries no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
