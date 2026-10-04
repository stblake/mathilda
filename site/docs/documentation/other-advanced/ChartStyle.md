# ChartStyle

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ChartStyle`**

BarChart/Histogram option: color or list of colors cycling through bars. Defaults to the standard palette.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= ChartStyle -> {Red, Blue}
Out[1]= ChartStyle -> {RGBColor[1, 0, 0], RGBColor[0, 0, 1]}
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[ChartStyle], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[BarChart[{3, 1, 2}, ChartStyle -> {Red, Green, Blue}]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `ChartStyle` is an inert option keyword for `BarChart` and `Histogram`: it gives a colour, or a list of colours cycled through the bars; it defaults to the standard palette. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `ChartStyle -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `ChartStyle -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `BarChart` and `Histogram`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `ChartStyle` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ChartStyle` is an inert option keyword for `BarChart` and `Histogram` — a colour or list of colours cycled through the bars. It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ChartStyle -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
