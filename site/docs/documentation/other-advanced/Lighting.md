# Lighting

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Lighting`**

A Graphics3D / Plot3D / ParametricPlot3D option controlling surface shading. Lighting -\> Automatic (default): per-face Lambertian (flat) shading with a fixed directional light; ambient 0.3, diffuse 0.7. Lighting -\> None (or False): disables shading and draws surfaces in their raw PlotStyle/ColorFunction color.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= Lighting -> None
Out[1]= Lighting -> None
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[Lighting], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, Lighting -> None]]
Out[3]= Graphics3D
```

## Implementation notes

**Definition.** `Lighting` is an inert option keyword for `Graphics3D`, `Plot3D` and `ParametricPlot3D`: it controls surface shading for 3D graphics; Automatic (the default) is per-face Lambertian shading with a fixed directional light, while None (or False) disables shading and draws the raw PlotStyle/ColorFunction colour. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `Lighting -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `Lighting -> value` rule by the owning plot
builtin while it assembles the `Graphics3D[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `Graphics3D`, `Plot3D` and `ParametricPlot3D`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `Lighting` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Lighting` is an inert option keyword for `Graphics3D`, `Plot3D` and `ParametricPlot3D` — surface shading for 3D graphics (Automatic = per-face Lambertian, None/False = flat raw colour). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `Lighting -> value` option rule while the owning builtin builds its `Graphics3D[...]`
result. The right-hand value evaluates normally before the builtin sees it.
