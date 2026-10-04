# Red

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Red`**

The named colour RGBColor\[1, 0, 0\].

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

The named primary evaluates straight to its RGBColor literal

```mathematica
In[1]:= Red
Out[1]= RGBColor[1, 0, 0]
```

An OwnValue, not an inert head -- it has already resolved

```mathematica
In[2]:= FullForm[Red]
Out[2]= RGBColor[1, 0, 0]
```

Usable wherever a colour directive is expected

```mathematica
In[3]:= Head[Graphics[{Red, Disk[]}]]
Out[3]= Graphics
```

An option value is a resolved colour too

```mathematica
In[4]:= FullForm[FrameStyle -> Red]
Out[4]= Rule[FrameStyle, RGBColor[1, 0, 0]]
```

## Implementation notes

**Definition.** `Red` is one of Mathilda's named colour constants. It is not an inert
head but an **OwnValue**: `register_color("Red", rgb_color(1, 0, 0))` in `graphics_init`
(`src/graphics/graphics_init.c`) binds the bare symbol `Red` to the literal
`RGBColor[1, 0, 0]` and stamps it `Protected`. The docstring is kept centrally in
`src/info.c`, beside `CMYKColor`, because `info_init` runs before `graphics_init`.

**Representation.** Evaluating `Red` yields `RGBColor[1, 0, 0]`, so `FullForm[Red]` is
`RGBColor[1, 0, 0]` and `Head[Red]` is `RGBColor` — the OwnValue has already fired by the
time anything inspects the result. The helper `color_component` prints an exact `0` or `1`
as a plain integer (matching Mathematica's own `InputForm[Red]`) and any other value as a
real. The saturated primaries use `RGBColor`; the greyscale names (`Black`, `White`,
`Gray`, `LightGray`) use `GrayLevel`.

**Usage & limits.** Because it resolves to a genuine colour literal wherever a bare
argument is evaluated, `Red` works anywhere a style directive is expected —
`Graphics[{Red, Disk[]}]`, `PlotStyle -> Red`, an `ArrayPlot` cell, `Directive[Red, Thick]`
— exactly as the long `RGBColor[...]` form would. `Plot` evaluates each option value once
before storing it (`split_options` in `plot.c`), so `Epilog -> {Red, ...}` also sees the
resolved literal. `Red` carries no opacity; use `Opacity[a]` or `RGBColor[r, g, b, a]` for
a translucent colour.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Red` is the named colour `RGBColor[1, 0, 0]`. It is bound as a `Protected` OwnValue, not
as an inert head, so it resolves to a real colour literal wherever a bare argument is
evaluated: `FullForm[Red]` is `RGBColor[1, 0, 0]` and `Head[Red]` is `RGBColor`.

That is what lets it be dropped anywhere a style directive is accepted — inside
`Graphics`, as a `PlotStyle`/`Epilog`/`FrameStyle` value, or in a `Directive`. The three
exact-`0`/`1` components print as integers, matching Mathematica's `InputForm`. `Red`
carries no opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]` for a
translucent colour.
