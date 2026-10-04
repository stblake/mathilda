---
references:
  - "R. M. Haralick and L. G. Shapiro, *Computer and Robot Vision*, vol. 1 (Addison-Wesley, 1992), ch. 2 (connected components labelling)."
source: src/imagefilter.c
---
**Algorithm.** `builtin_morphologicalcomponents` labels the connected components
of the foreground by **two-pass union-find**, returning an integer matrix with
background `0` and components numbered `1..k` in raster order of first
appearance. `MorphologicalComponents[image, t]` takes pixels strictly above `t`
as foreground (default `0`); `CornerNeighbors -> False` switches from the default
8-connectivity to 4. The **first pass** walks in raster order and can only see
already-visited neighbours (W, NW, N, NE for 8-connectivity), assigning
provisional labels and recording equivalences with `cc_union`/`cc_find` (path
compression; union by lower index, kept deterministic rather than by rank since
the second pass relabels anyway). A U-shaped region is why one pass is
insufficient — its two arms get different labels that the base reveals equal. The
**second pass** resolves each pixel to its representative and **relabels to
`1..k`** in raster order of first appearance, so there are no gaps and
`Max[result]` is exactly the component count. Connectivity is the one
discriminating property: two pixels touching only at a corner are **one**
component under 8 and **two** under 4.

**Data structures.** A `size_t` parent array for union-find over the pixel grid
and a label matrix; the result is a plain nested-`List` **integer matrix**, *not*
an `Image` — deliberately, because `Image` type inference would call a label
array of `1..12` a `"Byte"` image and `ImageData` would then divide every label
by 255. Labels are indices, not brightnesses, and must not be scaled.

**Complexity / limits.** `O(width · height · α(width · height))` — effectively
linear in the pixel count with path-compressed union-find. Foreground is
`pixel > t`; a colour image's pixels are compared on their stored channel values.
