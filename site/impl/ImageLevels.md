---
source: src/imagefilter.c
---
**Algorithm.** `builtin_imagelevels` returns the pixel histogram as **data**, a
list of `{level, count}` pairs — not a plot (`Histogram` over the result gives
the picture). `ImageLevels[image]` chooses the bin count from the type: a `"Bit"`
image uses its 2 natural levels, a `"Byte"` image its 256, because those *are*
the distinct values; only a `"Real"` image, which has no natural set, is binned
into 256 over `[0, 1]` (`IMG_LEVEL_BINS`). `ImageLevels[image, n]` forces `n`
bins. It loads the whole buffer (plane or volume), clamps each value to `[0, 1]`,
maps it to a bin with `floor(v · (nbins − 1) + 0.5)`, and accumulates a count.
The returned level for bin `i` is `i/(nbins − 1)`, on the **same unit scale**
`ImageData` uses, so a level is directly comparable to a pixel value without
rescaling. Every pixel lands in exactly one bin, so the counts sum to the pixel
count exactly.

**Data structures.** A `calloc`'d `double` count array of `nbins`, filled by one
pass over the flat unit-scale buffer (`image_load` / `image3d_load`), then
marshalled into a `List` of two-element `{real level, integer count}` pairs.

**Complexity / limits.** `O(n + nbins)`, `n` the pixel-channel count. `n` bins
must be between `2` and `65536`. Accepts volumes as well as planes; a colour
image's channels are binned together into the same histogram.
