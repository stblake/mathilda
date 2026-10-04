# Cyan

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Cyan`**

The named colour RGBColor\[0, 1, 1\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The named secondary evaluates straight to its RGBColor literal

```mathematica
In[1]:= Cyan
Out[1]= RGBColor[0, 1, 1]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[Cyan]
Out[2]= RGBColor[0, 1, 1]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{Cyan, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `Cyan` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("Cyan", rgb_color(0, 1, 1))` in `graphics_init`
(`src/graphics/graphics_init.c`) makes the bare symbol `Cyan` evaluate to the literal
`RGBColor[0, 1, 1]` and marks it `Protected`. The docstring is kept centrally in
`src/info.c`, beside `CMYKColor`.

**Representation.** `Cyan` evaluates to `RGBColor[0, 1, 1]`, so `FullForm[Cyan]` is
`RGBColor[0, 1, 1]` and `Head[Cyan]` is `RGBColor`. The `color_component` helper prints the
exact `0`/`1` components as integers. `Cyan`, `Magenta` (`RGBColor[1, 0, 1]`) and `Yellow`
(`RGBColor[1, 1, 0]`) are the secondary (subtractive) primaries — the same colour a
`CMYKColor` single-ink directive converts to.

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Cyan, Disk[]}]`, `PlotStyle -> Cyan`, `Directive[Cyan, Thick]`. `Plot` evaluates
its option values once before storing them, so `Epilog -> Cyan` also sees
`RGBColor[0, 1, 1]`. It carries no opacity; use `Opacity[a]` or the four-argument
`RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Cyan` is the named colour `RGBColor[0, 1, 1]`. It is bound as a `Protected` OwnValue, so
it resolves to a real colour literal wherever a bare argument is evaluated: `FullForm[Cyan]`
is `RGBColor[0, 1, 1]` and `Head[Cyan]` is `RGBColor`.

`Cyan`, `Magenta` (`RGBColor[1, 0, 1]`) and `Yellow` (`RGBColor[1, 1, 0]`) are the
secondary (subtractive) primaries — the colours a single-ink `CMYKColor` directive maps
to. The exact-`0`/`1` components print as integers. `Cyan` carries no opacity; use
`Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
