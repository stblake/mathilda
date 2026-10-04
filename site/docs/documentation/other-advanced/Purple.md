# Purple

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Purple`**

The named colour RGBColor\[0.5, 0, 0.5\].

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A named colour evaluates to its colour directive

```mathematica
In[1]:= Purple
Out[1]= RGBColor[0.5, 0, 0.5]
```

Same thing, confirmed

```mathematica
In[2]:= Purple === RGBColor[0.5, 0, 0.5]
Out[2]= True
```

Its components

```mathematica
In[3]:= List @@ Purple
Out[3]= {0.5, 0, 0.5}
```

## Implementation notes

**Definition.** `Purple` is one of Mathilda's **named colour constants** — purple. Like
Mathematica's own named colours it is not an inert head but a symbol whose OwnValue is
an `RGBColor` literal, so it resolves to a real colour directive anywhere one is expected:
`Purple` evaluates to `RGBColor[0.5, 0, 0.5]`.

**Representation.** `register_color("Purple", ...)` in `src/graphics/graphics_init.c`
installs the OwnValue `Purple -> RGBColor[0.5, 0, 0.5]` and marks the symbol `Protected`; the
component builder prints an exact `0` or `1` as a plain integer and anything else as
a real, matching Mathematica's own `InputForm`. The docstring lives centrally in
`src/info.c`. Because it is an OwnValue (not an inert head), a non-held graphics
argument or an evaluated-once plot option sees the expanded `RGBColor[0.5, 0, 0.5]`, exactly as a
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

`Purple` is one of Mathilda's named colour constants. Rather than being an inert head, it
is a symbol whose OwnValue is the directive `RGBColor[0.5, 0, 0.5]`, so it resolves to a real colour
wherever one is expected — as a `PlotStyle` / `VectorStyle` directive, inside
`Graphics` primitives, or in `ColorRules`.

It expands to an `RGBColor` directive, like most of the named colours. The exact components `0` and `1` print as plain integers (matching
Mathematica's `InputForm`). `Purple` is `Protected`.
