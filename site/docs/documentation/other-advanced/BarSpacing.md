# BarSpacing

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`BarSpacing`**

BarChart/Histogram option: gap between bars as a fraction of bar width (default 0.2). 0 = touching bars; 1 = all gap.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

An inert option keyword: it just names the left of a rule

```mathematica
In[1]:= BarSpacing -> 0.3
Out[1]= BarSpacing -> 0.3
```

Protected, with no value of its own

```mathematica
In[2]:= MemberQ[Attributes[BarSpacing], Protected]
Out[2]= True
```

Accepted by its owning plot builtin

```mathematica
In[3]:= Head[BarChart[{3, 1, 2}, BarSpacing -> 0.5]]
Out[3]= Graphics
```

## Implementation notes

**Definition.** `BarSpacing` is an inert option keyword for `BarChart` and `Histogram`: it sets the gap between bars as a fraction of bar width (default 0.2; 0 gives touching bars, 1 all gap). It has
**no builtin and no OwnValue** — only the `Protected` attribute bit is set on the
symbol in `graphics_init.c`, where its docstring also lives. It is a name, not a
value, and exists so that `BarSpacing -> ...` can be written in an option list.

**Representation.** A bare `Protected` `EXPR_SYMBOL`. It never evaluates on its own;
it is read as the left-hand side of a `BarSpacing -> value` rule by the owning plot
builtin while it assembles the `Graphics[...]` result, through the shared options
plumbing (`OptionValue`-style lookup). The value on the right is an ordinary
expression and is evaluated as usual — so e.g. a colour name on the right reaches
the builtin already reduced to its `RGBColor[...]` literal.

**Usage & limits.** Meaningful only inside a call to `BarChart` and `Histogram`; elsewhere it is an
unused symbol. An unrecognised value is handled (or defaulted) by the owning
builtin, not by `BarSpacing` itself, which carries no logic of its own.

**Attributes:** `Protected`.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`BarSpacing` is an inert option keyword for `BarChart` and `Histogram` — the gap between bars as a fraction of bar width (default 0.2; 0 gives touching bars, 1 all gap). It has no builtin and no
own-value; it is `Protected` and stays symbolic, read only as the left-hand side of
a `BarSpacing -> value` option rule while the owning builtin builds its `Graphics[...]`
result. The right-hand value evaluates normally before the builtin sees it.
