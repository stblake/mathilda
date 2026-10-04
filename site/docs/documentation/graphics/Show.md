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

## Examples (8)

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

### Applications (4)

```mathematica
In[5]:= Show[Graphics[{Red, Line[{{0, 0}, {1, 1}}]}], Graphics[Point[{2, 2}]]][[1]]
Out[5]= {{RGBColor[1, 0, 0], Line[{{0, 0}, {1, 1}}]}, {Point[{2, 2}]}}

In[6]:= Head[Show[Plot[Sin[x], {x, 0, 1}], Plot[Cos[x], {x, 0, 1}]]]
Out[6]= Graphics

In[7]:= Cases[Show[Graphics[Line[{{0, 0}, {1, 1}}], PlotRange -> {{0, 1}, {0, 1}}], Graphics[Point[{2, 3}]]], (PlotRange -> v_) :> v]
Out[7]= {{{0.0, 2.0}, {0.0, 3.0}}}

In[8]:= Show[Graphics[{Point[{0, 0}]}], Axes -> True]
Out[8]= -Graphics-
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

**Algorithm.** `builtin_show` is the engine's compositor and the REPL's display
entry point: any top-level result whose head is `Graphics`/`Graphics3D` is
auto-rendered, and `Show` is what normalises and merges those objects. It is not
`HoldAll` — its inputs are already-evaluated graphics. It partitions arguments
into leading graphics (descending into `List`s via `collect_graphics`) and
trailing options, and **declines** (returns `NULL`) if there are no graphics, if
a graphics follows an option, or if an argument is neither a graphics nor a list
of them; 2-D and 3-D objects cannot mix (all inputs must share one container
head). `Show[g]` / `Show[g, opts]` (`show_single`) keeps `g`'s container head and
merges each option by name, later-wins, appending unknown ones. `Show[g1, g2,
...]` (`show_multi`) overlays: each input's primitives become **one `List`
directive scope** (`scoped_prims`) so a colour or `Dashing` in one input cannot
restyle another, and for 2-D an input's own `PlotStyle` (first entry of a list)
is baked into its scope. Options are taken from `g1` (except `PlotStyle`,
`PlotRange`, and the stale `$PlotResample`/`$PlotLegendData` nodes) unless the
`Show` options override them. `PlotRange` is the **union** (`merged_plot_range`):
each input contributes an explicit numeric range if it has one, else the extent
of its primitives (`prims_extent` over `Line`/`Point`/`Polygon`/`Arrow`/
`Rectangle`/`Disk`/`Circle`/`Text`), and the union is emitted only if some input
fixed a range; otherwise it stays `Automatic`. `$PlotLegendData` nodes are
concatenated; `$PlotResample` is dropped from a combination.

**Data structures.** `PtrVec` (a borrowed-pointer vector) for the collected
graphics and option lists; `Box {lo[3], hi[3], any}` with `SHOW_MAX_DIM = 3` for
the range union. Under a `USE_GRAPHICS=0` build `show.c` carries a text-placeholder
stub that still performs the option merge and returns the merged object.

**Complexity / limits.** Linear in the primitive count for the range scan.
`prims_extent` reads only a fixed whitelist of primitive heads, so a head outside
it contributes no extent; a combined object keeps only `g1`'s non-merged options
(the other inputs' options, apart from `PlotStyle`/legends, are discarded).
`$RaylibVerbose` gates the backend's trace log.

- `Protected`.
- Declines to evaluate (stays unevaluated) if its first argument isn't a
  `Graphics[...]` expression, or any trailing argument isn't a `Rule`.
- When Raylib isn't compiled in, prints a one-line message instead of
  opening a window; the option merge still happens and the merged
  `Graphics[...]` is still returned.

**Attributes:** `Protected`.

## References

**See also:** [$RaylibVerbose](../../other-advanced/$RaylibVerbose/), [Plot](../../graphics/Plot/), [ListPlot](../../graphics/ListPlot/), [ParametricPlot](../../graphics/ParametricPlot/), [Plot3D](../../graphics/Plot3D/), [Animate](../../graphics/Animate/), [Manipulate](../../graphics/Manipulate/), [E](../../mathematical-constants/E/)

- Source: [`src/graphics/show.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/show.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

### Notes

`Show` is the engine's compositor and, via the REPL front end, its display path:
any top-level `Graphics`/`Graphics3D` result is auto-rendered, so `Plot[...]`,
`Show[...]` and `g // Graphics` all reach the window through one route. `Show[g,
opts]` merges options into `g` (later-wins, unknowns appended); `Show[g1, g2,
...]` overlays, wrapping **each input's primitives in its own `List` directive
scope** so a colour or `Dashing` in one input cannot restyle another, and baking
each 2-D input's own `PlotStyle` into its scope.

`PlotRange` of a combination is the **union** of the inputs' ranges — an input
without an explicit range contributes the extent of its primitives — and stays
`Automatic` unless at least one input fixed a range. Options other than the merged
ones are taken from the first graphic. `Show` declines (stays unevaluated) if its
first argument is not a graphics, if a graphics follows an option, or if 2-D and
3-D graphics are mixed.
