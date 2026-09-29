# Show

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Show[graphics, opts...]`**

Displays graphics (a Graphics\[...\] object) and returns it, with the options opts overriding its own.

**`Show[g1, g2, ..., opts...]`**

Combines several graphics (Plot, ListPlot, Graphics, ... outputs) into one: their primitives are overlaid, each in its own style scope; options come from g1 unless given in opts; the plot range is the union of the inputs' ranges.

**`Show[{g1, g2, ...}, opts...]`**

The same, for a list of graphics.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Show[Graphics[{Red, Line[{{0, 0}, {1, 1}}]}], Graphics[Point[{2, 2}]]][[1]]
Out[1]= {{RGBColor[1, 0, 0], Line[{{0, 0}, {1, 1}}]}, {Point[{2, 2}]}}

In[2]:= Head[Show[Plot[Sin[x], {x, 0, 1}], Plot[Cos[x], {x, 0, 1}]]]
Out[2]= Graphics
```

### Options (2)

```mathematica
In[3]:= Cases[Show[Graphics[Line[{{0, 0}, {1, 1}}], PlotRange -> {{0, 1}, {0, 1}}], Graphics[Point[{2, 3}]]], (PlotRange -> v_) :> v]
Out[3]= {{{0.0, 2.0}, {0.0, 3.0}}}

In[4]:= Show[Graphics[{Point[{0,0}]}], Axes -> True]
Out[4]= -Graphics-
```

## Algorithm

show.c — Show[graphics, opts...] / Show[g1, g2, ..., opts...] plus the no-Raylib graphics_show() stub. When USE_GRAPHICS is compiled in, render.c provides the real graphics_show(); this file's stub is excluded then (see the #ifndef below) so there's exactly one definition either way.

Show follows Mathematica:

```text
  Show[g, opts]           g with opts overriding its own options.
  Show[g1, g2, ..., opts] the graphics overlaid. Each input's primitives
                          become one List -- a directive scope -- so a
                          colour or Dashing in g1 never restyles g2. An
                          input's own PlotStyle (the style a single-curve
                          Plot/ListPlot is drawn in) is baked into its
                          scope, because the combined object keeps only
                          g1's options. Options come from g1 unless opts
                          overrides them; PlotRange is the union of the
                          inputs' ranges (a graphic without an explicit
                          range contributes its primitives' extent), and
                          is left automatic when no input fixes one.
  Show[{g1, g2, ...}, ..] the same for (nested) lists of graphics.
```

Plot's $PlotResample metadata re-samples only its own curves, so it is dropped from a combination (the static primitives stay); the inputs' $PlotLegendData entries are concatenated into one legend. 2D and 3D graphics cannot be mixed (the call stays unevaluated).

## Implementation notes

- `Protected`.
- Declines to evaluate (stays unevaluated) if its first argument isn't a
  `Graphics[...]` expression, or any trailing argument isn't a `Rule`.
- When Raylib isn't compiled in, prints a one-line message instead of
  opening a window; the option merge still happens and the merged
  `Graphics[...]` is still returned.

**Attributes:** `Protected`.

## References

**See also:** [$RaylibVerbose](../../other-advanced/$RaylibVerbose/), [Plot](../../graphics/Plot/), [ListPlot](../../graphics/ListPlot/), [ParametricPlot](../../graphics/ParametricPlot/), [Plot3D](../../graphics/Plot3D/), [Animate](../../graphics/Animate/), [Manipulate](../../graphics/Manipulate/), [E](../../mathematical-constants/E/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
