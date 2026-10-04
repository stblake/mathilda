---
references:
  - "T. Y. Zhang and C. Y. Suen, *A fast parallel algorithm for thinning digital patterns*, Comm. ACM **27** (1984) 236-239."
source: src/imagethin.c
---
**Algorithm.** `builtin_pruning` removes one pixel from every **free end** of the
foreground. It builds a foreground mask (`mask_from_image`, channels averaged and
thresholded at `0.5`) and runs `prune_pass` `n` times — once by default, `n`
times for `Pruning[image, n]`, which shortens each branch by up to `n` and
deletes any branch shorter than that; `Pruning[image, 0]` is the image unchanged.
Each pass marks an **end point** — a foreground pixel with exactly **one**
foreground neighbour — and deletes all marked pixels together after the pass
(mirroring `Thinning`, so one deletion does not change a neighbour's verdict
mid-pass). The "exactly one neighbour" rule is deliberate: an **isolated** pixel
has *zero* foreground neighbours and so is not an end point and survives — pruning
shortens branches, it does not erase specks, and a rule that deleted isolated
pixels would quietly remove every one-pixel component. Pruning is the companion
to `Thinning`, used to clean the short spurs a skeleton grows at boundary
irregularities. The result is a `"Bit"` image.

**Data structures.** An `unsigned char` foreground mask of `width · height` and a
per-pass `calloc`'d deletion map. The 8-neighbourhood count reuses the same
`neighbours`/counting helpers as `Thinning`. Built into a `"Bit"` image
(`image_build_bit`).

**Complexity / limits.** `O(n · width · height)` for `n` passes. Operates on a
binary plane; a non-binary image is thresholded at `0.5`.
