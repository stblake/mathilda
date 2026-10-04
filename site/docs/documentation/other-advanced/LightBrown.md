# LightBrown

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightBrown`**

The named colour RGBColor\[0.94, 0.91, 0.88\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A pale pastel -- its own constant, not a tint of Brown

```mathematica
In[1]:= LightBrown
Out[1]= RGBColor[0.94, 0.91, 0.88]
```

An OwnValue, not an inert head -- already resolved

```mathematica
In[2]:= FullForm[LightBrown]
Out[2]= RGBColor[0.94, 0.91, 0.88]
```

Usable wherever a directive is expected

```mathematica
In[3]:= Head[Graphics[{LightBrown, Disk[]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `LightBrown` is one of Mathilda's named colour constants, bound as an
**OwnValue**: `register_color("LightBrown", rgb_color(0.94, 0.91, 0.88))` in
`graphics_init` (`src/graphics/graphics_init.c`) makes the bare symbol evaluate to the
literal `RGBColor[0.94, 0.91, 0.88]` and marks it `Protected`. The docstring is kept
centrally in `src/info.c`.

**Representation.** `LightBrown` evaluates to `RGBColor[0.94, 0.91, 0.88]`, so
`FullForm[LightBrown]` is `RGBColor[0.94, 0.91, 0.88]` and `Head[LightBrown]` is
`RGBColor`. Its components are machine reals (none is exactly `0` or `1`). It is a
*separate* constant from `Brown`, not a computed tint of it — the `Light*` family is its
own set of pale pastel entries (`LightRed`, `LightGreen`, `LightBlue`, …).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{LightBrown, Disk[]}]`, `PlotStyle -> LightBrown`,
`Directive[LightBrown, Thick]` — and as a pale background-style fill. It carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightBrown` is the named colour `RGBColor[0.94, 0.91, 0.88]`. It is bound as a `Protected`
OwnValue, so it resolves to a real colour literal wherever a bare argument is evaluated:
`FullForm[LightBrown]` is `RGBColor[0.94, 0.91, 0.88]` and `Head[LightBrown]` is `RGBColor`.

It is a *separate* constant from `Brown`, not a computed tint — the `Light*` family is its
own set of pale pastel entries (`LightRed`, `LightGreen`, `LightBlue`, …). It carries no
opacity; use `Opacity[a]` or the four-argument `RGBColor[r, g, b, a]`.
