# Green

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Green`**

The named colour RGBColor\[0, 1, 0\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The named primary evaluates straight to its RGBColor literal

```mathematica
In[1]:= Green
Out[1]= RGBColor[0, 1, 0]
```

An OwnValue, not an inert head -- it has already resolved

```mathematica
In[2]:= FullForm[Green]
Out[2]= RGBColor[0, 1, 0]
```

Usable wherever a directive is

```mathematica
In[3]:= Head[Graphics[{Green, Line[{{0, 0}, {1, 1}}]}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `Green` is one of Mathilda's named colour constants, bound as an
**OwnValue** rather than registered as an inert head: `register_color("Green",
rgb_color(0, 1, 0))` in `graphics_init` (`src/graphics/graphics_init.c`) makes the bare
symbol `Green` evaluate to the literal `RGBColor[0, 1, 0]` and marks it `Protected`. The
docstring lives centrally in `src/info.c` (set in `info_init`, which runs first).

**Representation.** `Green` evaluates to `RGBColor[0, 1, 0]`, so `FullForm[Green]` is
`RGBColor[0, 1, 0]` and `Head[Green]` is `RGBColor`. The `color_component` helper emits an
exact `0`/`1` as an integer and anything else as a real, matching Mathematica's
`InputForm`. `Green` is one of the saturated primaries built with `RGBColor` (the
greyscale names use `GrayLevel`).

**Usage & limits.** As a resolved colour literal it is accepted anywhere a directive is —
`Graphics[{Green, Line[...]}]`, `PlotStyle -> Green`, `Directive[Green, Dashed]`. `Plot`
evaluates its option values once before storing them, so a held option like
`Epilog -> Green` also sees `RGBColor[0, 1, 0]`. It carries no opacity; use `Opacity[a]` or
the four-argument `RGBColor` for translucency.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Green` is the named colour `RGBColor[0, 1, 0]`. It is bound as a `Protected` OwnValue, so
it resolves to a real colour literal wherever a bare argument is evaluated: `FullForm[Green]`
is `RGBColor[0, 1, 0]` and `Head[Green]` is `RGBColor`.

That resolution is what lets it be used anywhere a style directive is accepted — inside
`Graphics`, as a `PlotStyle`/`Epilog` value, or in a `Directive`. The exact-`0`/`1`
components print as integers. `Green` carries no opacity; use `Opacity[a]` or the
four-argument `RGBColor[r, g, b, a]`.
