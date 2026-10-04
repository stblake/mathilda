---
references:
  - "P. Heckbert, *Color image quantization for frame buffer display*, Computer Graphics (SIGGRAPH) **16** (1982) 297-307."
source: src/imagecolor.c
---
**Algorithm.** `builtin_colorquantize` reduces an image to at most `n` colours by
**median cut**. It loads the pixels to a flat buffer, copies each as an RGB
triple into `qbuf` (a grey pixel is replicated across the three channels), and
starts with one box covering all pixel indices. While there are fewer than `n`
boxes it repeatedly splits the box with the **widest single-channel spread** —
widest spread rather than most pixels, since a large box of nearly identical
colours does not need splitting and a small one spanning half the spectrum does.
To split, it `qsort`s that box's indices on the widest channel (ties broken by
index, so the palette does not depend on qsort's internal choices) and cuts at
the median. When every box is a single colour the loop stops early. Finally each
box collapses to its **mean** colour, written back to every pixel it owns; a grey
output takes the mean's luminance, and an alpha channel passes through unchanged.
Median cut is chosen over k-means precisely because it is **deterministic** — a
palette depending on a random seed could not be tested or documented.

**Data structures.** `qbuf` is `n_px × 3` doubles (the colours), `qidx` an index
permutation sorted per box, `boxes` an array of `{start, len}` ranges into
`qidx`, and `out` the `n_px × channels` result buffer. File-scope `qbuf`/`qidx`/
`qchan` are the comparator's context for `qsort`. The result is a `"Real"` image.

**Complexity / limits.** Each of the up-to-`n` splits sorts a box:
`O(n · n_px log n_px)` worst case, with the per-split scan for the widest channel
`O(n · n_px)`. `n` is `1..4096`. The channel count is preserved and alpha passes
through.
