---
source: src/imagecompose.c
---
**Algorithm.** `builtin_imagecompose` alpha-composites `over` onto `base`. The result keeps
the base's size — composition is "draw on this", not "make something bigger" — and any part of
`over` falling outside is clipped. The overlay is centred by default; `ImageCompose[base, over,
{x, y}]` places its centre at image coordinate `{x, y}` (x from the left, **y from the
bottom**), which the storage-order flip `oy = base.h - py - over.h/2` converts to a row offset.
`ImageCompose[base, {over, opacity}]` scales the overlay's alpha by a constant.

Per destination pixel the standard over operator runs: `out = o·α + b·(1−α)` for each colour
channel, with `α` the overlay's own alpha at that point times the opacity. Channel counts are
reconciled once by `promote`, which combines in the larger colour count — a grey overlay on a
colour base produces colour, since grey means "the same in every channel" and is replicated
rather than zero-padded. The result carries an alpha channel only if the base did
(`out_alpha = base_alpha + α·(1 − base_alpha)`); compositing onto an opaque image stays opaque.

**Data structures.** Two decoded `Img` buffers; a fresh
`base.w · base.h · out_channels` buffer, wrapped by `image_build_real` as a packed `"Real"`
image (every image head returns a packed buffer — `make check-image-packing`). `sample`
handles grey-to-colour replication and absent alpha in one accessor.

**Complexity / limits.** `O(base.w · base.h · channels)`. Placement outside the base is
clipped, not an error; the opacity constant must lie in `[0, 1]`.
