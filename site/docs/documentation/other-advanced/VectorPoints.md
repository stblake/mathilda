# VectorPoints

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`VectorPoints`**

VectorPlot option: integer n specifies an n×n seed grid (default 15). Automatic also uses 15.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= VectorPoints -> 20
Out[1]= VectorPoints -> 20
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[VectorPoints], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 5]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `VectorPoints` is an inert option keyword for `VectorPlot`: it gives the integer n specifying an n*n seed grid of arrows (default 15; Automatic also uses 15). It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `VectorPoints -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `VectorPoints -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `VectorPlot`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `VectorPoints` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`VectorPoints` is an inert option keyword for `VectorPlot` — an integer n giving an n*n seed grid of arrows (default 15). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `VectorPoints -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
