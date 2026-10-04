---
source: src/graphics/show.c
---
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
