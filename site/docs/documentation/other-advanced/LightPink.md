# LightPink

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightPink`**

The named colour RGBColor\[1, 0.925, 0.925\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A pale pastel -- distinct from Pink

```mathematica
In[1]:= LightPink
Out[1]= RGBColor[1, 0.925, 0.925]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[LightPink]
Out[2]= RGBColor[1, 0.925, 0.925]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{LightPink, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `LightPink` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightPink", rgb_color(1, 0.925, 0.925))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the literal
`RGBColor[1, 0.925, 0.925]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `LightPink` evaluates to `RGBColor[1, 0.925, 0.925]`, so
`FullForm[LightPink]` is `RGBColor[1, 0.925, 0.925]` and `Head[LightPink]` is `RGBColor`.
The exact `1` red component prints as an integer and the `0.925` green/blue as reals. It
is a distinct constant from `Pink` (`RGBColor[1, 0.5, 0.5]`), one of the pale `Light*`
pastel family.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightPink, Disk[]}]`, `PlotStyle -> LightPink`, or a pale fill/background. It
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightPink` is the named colour `RGBColor[1, 0.925, 0.925]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightPink]` is `RGBColor[1, 0.925, 0.925]` and `Head[LightPink]` is `RGBColor`.

The exact `1` red component prints as an integer and the `0.925` green/blue as reals. It
is a distinct constant from `Pink` (`RGBColor[1, 0.5, 0.5]`), one of the pale `Light*`
pastel family. It carries no opacity; use `Opacity[a]` or the four-argument
`RGBColor[r, g, b, a]`.
