---
source: src/imagecompose.c
---
**Algorithm.** `builtin_imageassemble` tiles a grid of images into one.
`ImageAssemble[{{a, b}, {c, d}}]` is a 2-D grid; `ImageAssemble[{a, b}]` a single row (a flat
list is one row, the reading order a list already implies — distinguished by probing whether
the first element is itself an image). Each tile is decoded with `image_load`. Row heights and
column widths are taken as the maximum over the row/column (`row_h[i]`, `col_w[j]`), so tiles
of different sizes line up on a grid; the total canvas is the sum of those. Each tile is
blitted at its natural size into its cell — nothing is stretched, because stretching would
resample an image the caller did not ask to resize, and any leftover gap is left blank.

The combined colour-channel count is the maximum over all tiles (`promote`), with grey
replicated into colour via `sample`. An alpha channel is added to the sheet iff **any** tile
had one, and where a sheet has alpha its gaps are transparent rather than opaque black — the
one property a sprite sheet must keep.

**Data structures.** An array of decoded `Img` tiles plus parallel `row_h`/`col_w` arrays; a
single `calloc`'d `total_w · total_h · out_channels` canvas, wrapped by `image_build_real` as a
packed `"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(total pixels)` to clear the canvas plus `O(sum of tile pixels)` to
copy; a ragged grid is accepted (short rows simply leave cells blank). A non-image element
makes the whole call decline.
