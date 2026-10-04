# Black

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Black`**

The named colour GrayLevel\[0\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A greyscale name resolves to GrayLevel, not RGBColor

```mathematica
In[1]:= Black
Out[1]= GrayLevel[0]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[Black]
Out[2]= GrayLevel[0]
```

Usable wherever a directive is

```mathematica
In[3]:= Head[Graphics[{Black, Line[{{0, 0}, {1, 1}}]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `Black` is one of Mathilda's named colour constants. It is bound as an
**OwnValue** — `register_color("Black", gray_color(0))` in `graphics_init`
(`src/graphics/graphics_init.c`) — so the bare symbol `Black` evaluates to the literal
`GrayLevel[0]` and is marked `Protected`. Unlike the saturated primaries it uses
`GrayLevel`, not `RGBColor`, matching Mathematica's own `InputForm[Black]`. The docstring
is kept centrally in `src/info.c`.

**Representation.** `Black` evaluates to `GrayLevel[0]`, so `FullForm[Black]` is
`GrayLevel[0]` and `Head[Black]` is `GrayLevel`. The single component is an exact integer
`0` (the `color_component` helper prints an exact `0`/`1` as an integer). `Black`, `White`
(`GrayLevel[1]`), `Gray` (`GrayLevel[0.5]`) and `LightGray` (`GrayLevel[0.85]`) are the
greyscale family.

**Usage & limits.** As a resolved directive it is accepted wherever a colour is —
`Graphics[{Black, Line[...]}]`, `PlotStyle -> Black`, `FrameStyle -> Black`,
`Directive[Black, Thick]`. A renderer only ever has to understand the long `GrayLevel[...]`
form. It carries no opacity; use `Opacity[a]` or `GrayLevel[g, a]` for a translucent black.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Black` is the named colour `GrayLevel[0]`. Unlike the saturated primaries it uses
`GrayLevel`, not `RGBColor`, matching Mathematica's `InputForm`: `FullForm[Black]` is
`GrayLevel[0]` and `Head[Black]` is `GrayLevel`.

It is bound as a `Protected` OwnValue, so it resolves to a real colour literal wherever a
bare argument is evaluated, and can be dropped anywhere a style directive is accepted —
`Graphics`, `PlotStyle`, `FrameStyle`, `Directive`. Its companions are `White`
(`GrayLevel[1]`), `Gray` (`GrayLevel[0.5]`) and `LightGray` (`GrayLevel[0.85]`). It carries
no opacity; use `Opacity[a]` or `GrayLevel[g, a]`.
