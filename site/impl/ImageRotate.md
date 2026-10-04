---
source: src/imagegeom.c
---
**Algorithm.** `builtin_imagerotate` rotates counterclockwise (positive angle), defaulting to
a quarter turn. It accepts an angle in radians, an exact form such as `Pi/2` or `90 Degree`
(numericalised through `N` first, since `na_read_scalar` refuses a symbolic value), a side
(`Left`/`Right`/`Top`/`Bottom`, turning the top of the image to face it), or `side1 -> side2`.

The key split is between right angles and everything else. When the angle is a multiple of a
right angle (tested on the *angle*, `|q − round(q)| < 1e-9`), `rot90_run` takes a **pure index
permutation** — every pixel lands on another pixel's exact position, nothing is interpolated,
four quarter turns are exactly the identity, and an odd number of quarters swaps the
dimensions. Any other angle calls `rot_free_run`: a bilinear **inverse** map (iterate over
destination pixels, ask where each came from, so every output is filled exactly once —
forward mapping leaves holes where the rotation stretches). The inverse rotation uses
`[[c, −s], [s, c]]`, the transpose of the textbook y-up matrix, because image rows run *down*.
Area rotated in from outside the frame reads as 0, not the replicated edge, because that area
was never photographed and smearing the border across it would invent content.

**Data structures.** One decoded unit buffer in; a fresh buffer out (dimensions swapped for an
odd quarter turn), wrapped by `image_build_real` as a packed `"Real"` image (every image head
returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)`. A free angle keeps the input's dimensions (Mathematica's
`Full` would enlarge to enclose the rotated image), so corners rotate out of frame — an
intentional divergence.
