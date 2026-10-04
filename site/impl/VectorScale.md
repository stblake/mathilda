---
source: src/graphics/vectorplot.c
---
**Definition.** `VectorScale` is an **option symbol** for `VectorPlot` that controls how
arrow length encodes field magnitude. It has no builtin and no value of its own — it is
an inert keyword (`Protected`) read out of `VectorPlot`'s option sequence.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_VectorScale`). `VectorPlot`'s
option parser (`src/graphics/vectorplot.c`) evaluates the right-hand side of a
`VectorScale -> v` rule and sets an internal `scale_mode`: `Automatic` (the default) →
`VS_AUTOMATIC`, all arrows drawn at equal display length so only direction is shown;
`None` → `VS_NONE`, length proportional to field magnitude; a positive real `f` →
`VS_FRACTION` with arrow length `= f × grid spacing`. Any other value falls back to
`Automatic`.

**Usage & limits.** `Protected`. It carries meaning only inside a `VectorPlot` call; on
its own it just evaluates to itself. It has no entry in `Options[VectorPlot]`'s default
list — it is read directly from the given option sequence — so `SetOptions` on it has no
effect. Companion `VectorPlot` options are `VectorPoints` (seed-grid density) and
`VectorStyle` (global style directives).
