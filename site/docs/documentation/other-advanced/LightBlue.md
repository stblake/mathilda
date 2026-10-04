# LightBlue

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LightBlue`**

The named colour RGBColor\[0.87, 0.94, 1\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A named colour evaluates to its colour directive

```mathematica
In[1]:= LightBlue
Out[1]= RGBColor[0.87, 0.94, 1]
```

Same thing, confirmed

```mathematica
In[2]:= LightBlue === RGBColor[0.87, 0.94, 1]
Out[2]= True
```

Its components

```mathematica
In[3]:= List @@ LightBlue
Out[3]= {0.87, 0.94, 1}
```

## Implementation notes

**Definition.** `LightBlue` is one of Mathilda's **named colour constants** — a pale blue. Like
Mathematica's own named colours it is not an inert head but a symbol whose OwnValue is
an `RGBColor` literal, so it resolves to a real colour directive anywhere one is expected:
`LightBlue` evaluates to `RGBColor[0.87, 0.94, 1]`.

**Representation.** `register_color("LightBlue", ...)` in `src/graphics/graphics_init.c`
installs the OwnValue `LightBlue -> RGBColor[0.87, 0.94, 1]` and marks the symbol `Protected`; the
component builder prints an exact `0` or `1` as a plain integer and anything else as
a real, matching Mathematica's own `InputForm`. The docstring lives centrally in
`src/info.c`. Because it is an OwnValue (not an inert head), a non-held graphics
argument or an evaluated-once plot option sees the expanded `RGBColor[0.87, 0.94, 1]`, exactly as a
literal colour would.

**Usage & limits.** `Protected`. It is a value, not a function: feeding it a wrong
"argument" is meaningless. Use it as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`. The grey-scale constants
(`Black`, `White`, `Gray`, `LightGray`) expand to `GrayLevel`; the rest to
`RGBColor`.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`LightBlue` is one of Mathilda's named colour constants. Rather than being an inert head, it
is a symbol whose OwnValue is the directive `RGBColor[0.87, 0.94, 1]`, so it resolves to a real colour
wherever one is expected — as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`.

It expands to an `RGBColor` directive, like most of the named colours. The exact components `0` and `1` print as plain integers (matching
Mathematica's `InputForm`). `LightBlue` is `Protected`.
