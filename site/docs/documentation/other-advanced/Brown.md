# Brown

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Brown`**

The named colour RGBColor\[0.6, 0.4, 0.2\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A tertiary tone -- components are reals, not 0/1

```mathematica
In[1]:= Brown
Out[1]= RGBColor[0.6, 0.4, 0.2]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[Brown]
Out[2]= RGBColor[0.6, 0.4, 0.2]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{Brown, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `Brown` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("Brown", rgb_color(0.6, 0.4, 0.2))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol `Brown` evaluate to the literal
`RGBColor[0.6, 0.4, 0.2]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`.

**Representation.** `Brown` evaluates to `RGBColor[0.6, 0.4, 0.2]`, so `FullForm[Brown]`
is `RGBColor[0.6, 0.4, 0.2]` and `Head[Brown]` is `RGBColor`. Its components are not `0`
or `1`, so `color_component` renders each as a machine real — unlike the pure primaries,
whose components print as integers. `Brown` is one of the fuller tertiary tones (with
`Orange`, `Pink`, `Purple`).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Brown, Disk[]}]`, `PlotStyle -> Brown`, `Directive[Brown, Thick]`. `Plot`
evaluates its option values once before storing them, so `Epilog -> Brown` also sees the
resolved literal. The light companion `LightBrown` (`RGBColor[0.94, 0.91, 0.88]`) is a
separate constant, not a programmatic lightening of `Brown`. It carries no opacity; use
`Opacity[a]` or the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Brown` is the named colour `RGBColor[0.6, 0.4, 0.2]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[Brown]` is `RGBColor[0.6, 0.4, 0.2]` and `Head[Brown]` is `RGBColor`.

Its components are machine reals (none is exactly `0` or `1`), so they print as reals,
unlike the pure primaries. The light companion `LightBrown` (`RGBColor[0.94, 0.91, 0.88]`)
is a separate constant, not a programmatic lightening of `Brown`. `Brown` carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
