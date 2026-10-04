---
source: src/imagefilter.c
---
**Algorithm.** `builtin_erosion` is flat morphological erosion: the **minimum** over a
neighbourhood. `morph_builtin(res, MORPH_ERODE, false)` reads the structuring element with
`morph_element` — an integer radius `r` means a full `(2r+1)²` rectangle, a matrix means its
nonzero **support** (flat: the element's values do not enter the minimum, which keeps
`Erosion[img, BoxMatrix[1]]` and `Erosion[img, 1]` the same operation). Padding replicates the
edge, the same rule the convolutions use; this is what makes the morphology laws hold at the
border — zero padding would let an erosion see a black neighbour that is not there.

A full rectangle is **separable** for the minimum exactly as for a sum (the min over a
rectangle is the min over rows of the mins over columns), and a 1-D min is further reduced to
three comparisons per pixel independent of radius by **van Herk–Gil–Werman** (prefix/suffix
minima over blocks of width k), which is what makes erosion usable at large radii. The fast
path is bit-exact with the naive `morph_direct`, not approximate, because min is associative
and idempotent. Erosion is dual to Dilation: for a symmetric element `Erosion[f, k] =
1 − Dilation[1 − f, k]` exactly, which holds at the border only because replicate padding is
itself self-dual.

**Data structures.** A decoded unit buffer; the result through `image_build_typed`/`_real`. A
`"Bit"` image stays `"Bit"` (a min of stored values is a stored value), other types become
`"Real"`; either way it is a packed buffer (`make check-image-packing`).

**Complexity / limits.** `O(pixels)` for a full rectangle (van Herk), `O(pixels · |support|)`
for an arbitrary element. Radius ≤ 256 planar, ≤ 64 volumetric.
