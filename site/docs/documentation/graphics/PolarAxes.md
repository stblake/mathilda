# PolarAxes

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolarAxes`**

PolarPlot option: True requests a polar grid overlay (radial circles + angle labels). Currently accepted but not yet rendered; Cartesian axes are drawn instead.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= {Head[PolarAxes], MemberQ[Attributes[PolarAxes], Protected]}
Out[1]= {Symbol, True}
```

### Options (2)

```mathematica
In[2]:= PolarPlot[1, {t, 0, 2 Pi}, PolarAxes -> True]
Out[2]= -Graphics-

In[3]:= Cases[PolarPlot[Sin[2 t], {t, 0, 2 Pi}, PolarAxes -> True], (PolarAxes -> v_) :> v, Infinity]
Out[3]= {True}
```

### Applications (4)

```mathematica
In[4]:= PolarPlot[2, {t, 0, 2 Pi}, PolarAxes -> True]
Out[4]= -Graphics-

In[5]:= Cases[PolarPlot[Sin[2 t], {t, 0, 2 Pi}, PolarAxes -> True], (PolarAxes -> v_) :> v, Infinity]
Out[5]= {True}

In[6]:= MemberQ[Attributes[PolarAxes], Protected]
Out[6]= True

In[7]:= Head[PolarAxes]
Out[7]= Symbol
```

## Implementation notes

**Algorithm.** `PolarAxes` is not a plotter — it is an inert, `Protected` option
keyword for `PolarPlot`, registered in `graphics_init.c` with only an attribute
bit and a docstring (there is no builtin function and no OwnValue). A grep of the
whole tree finds `PolarAxes` referenced nowhere outside `graphics_init.c`: it is
not one of the names `split_options_param` consumes, so when passed as
`PolarAxes -> True` it falls into the generic option pass-through and rides
verbatim as a trailing `Rule` onto the resulting `Graphics[...]` object. The
renderer does not yet draw the requested polar grid (radial circles plus angular
degree/radian labels), so **Cartesian axes are drawn instead** — `PolarPlot`
injects the usual `Axes -> True`. The symbol therefore stays symbolic in every
context and has no evaluation rules of its own.

**Data structures.** None — it is a bare interned symbol carrying
`ATTR_PROTECTED` and a docstring; its value, when supplied, is stored as an
ordinary `Rule[PolarAxes, True]` option on the graphics object.

**Complexity / limits.** None. The feature is deliberately a documented
placeholder: the option is recognised and preserved so existing `PolarPlot`
calls that set it keep working, but the polar grid overlay is unimplemented and
the plot renders with Cartesian axes. This is recorded here rather than silently
dropped.

**Attributes:** `Protected`.

## References

**See also:** [PolarPlot](../../graphics/PolarPlot/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`PolarAxes` is not a plotter but an inert, `Protected` **option keyword** for
`PolarPlot`. `PolarAxes -> True` asks for a polar grid overlay — radial circles at
regular intervals plus angular degree/radian labels — but that overlay is not yet
rendered, so the plot is drawn with ordinary Cartesian axes instead. This is a
documented placeholder rather than a silent no-op: the option is recognised and
preserved (it passes through as a `Rule` onto the returned `Graphics[...]`), so
existing `PolarPlot` calls that set it keep working unchanged.

The symbol carries no builtin and no own-values; it has only the `Protected`
attribute and a docstring, and stays symbolic in every context.
