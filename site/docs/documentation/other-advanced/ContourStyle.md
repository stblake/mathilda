# ContourStyle

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ContourStyle`**

ContourPlot option: style directive(s) for contour lines. A list cycles through the levels; Automatic colours by height; None suppresses lines.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= ContourStyle -> Red
Out[1]= ContourStyle -> RGBColor[1, 0, 0]
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[ContourStyle], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[ContourPlot[x^2 + y^2, {x, -1, 1}, {y, -1, 1}, ContourStyle -> Red]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `ContourStyle` is an inert option keyword for `ContourPlot`: it gives the style directive(s) for the contour lines; a single directive applies to all levels, a list cycles through them, Automatic colours by height and None suppresses the lines. It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `ContourStyle -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `ContourStyle -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `ContourPlot`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `ContourStyle` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ContourStyle` is an inert option keyword for `ContourPlot` — the style directive(s) for contour lines; a list cycles through levels, Automatic colours by height, None suppresses lines. It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `ContourStyle -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
