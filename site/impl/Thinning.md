---
references:
  - "T. Y. Zhang and C. Y. Suen, *A fast parallel algorithm for thinning digital patterns*, Comm. ACM **27** (1984) 236-239."
source: src/imagethin.c
---
**Algorithm.** `builtin_thinning` reduces the foreground to a one-pixel-wide
skeleton by **Zhang–Suen** thinning. It builds a foreground mask
(`mask_from_image`: channels averaged, thresholded at `0.5`) and iterates
`thin_pass` until a pass deletes nothing (`Thinning[image]`) or for `n` passes
(`Thinning[image, n]`). Each full iteration runs **two subiterations**
(`thin_pass(..., second = false)` then `second = true`), and the two are what
preserve connectivity: deleting every individually-removable pixel in one pass
severs a diagonal line, since two diagonal neighbours can each be removable while
removing both disconnects the shape — the two conditions delete from opposite
sides on alternating passes. A pixel is marked for deletion when it has `2..6`
foreground neighbours (`neighbour_count`), exactly one `0→1` transition around
the 8-ring (`transitions`, so removing it would not break connectivity), and
satisfies the subiteration's pair of corner conditions. Crucially, marked pixels
are deleted **together after the whole pass** — deleting in place would let one
pixel's removal change the verdict on its neighbour mid-pass. Outside the image
counts as background. The result is always a `"Bit"` image and is always a subset
of the input.

**Data structures.** An `unsigned char` foreground mask of `width · height`, plus
a per-pass `unsigned char` deletion map (`calloc`) so marking and deletion are
separated. The eight neighbours are gathered in Zhang–Suen's order (`neighbours`,
P2 north then clockwise). Built into a `"Bit"` image (`image_build_bit`).

**Complexity / limits.** `O(width · height)` per subiteration; the number of
iterations to convergence is bounded by half the shape's thickness. A non-binary
image is thresholded at `0.5` — apply `Binarize` first for any other rule.
